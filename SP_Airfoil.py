# SOURCE PANEL METHOD - SINGLE AIRFOIL
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 01/01/19 - In MATLAB
# Updated   : 02/03/19 - Transferred from MATLAB to Python
#                      - Works as expected
#             02/09/20 - Added DAT airfoil loading option with XFOIL function
# Notes     : This code is not optimized, but is instead written in such a way
#             that it is easy to follow along with my YouTube video derivations
# 
# Functions Needed:
# - XFOIL.py
# - COMPUTE_IJ_SPM.py
# - STREAMLINE_SPM.py
# - COMPUTE_CIRCULATION.py
#
# Programs Needed:
# - xfoil.exe
# 
# Folder Needed:
# - Airfoil_DAT_Selig: folder containing all Selig-format airfoils
# 
# References
# - [1]: Panel Method Geometry
#           Link: https://www.youtube.com/watch?v=kIqxbd937PI
# - [2]: Normal Geometric Integral SPM, I(ij)
#           Link: https://www.youtube.com/watch?v=76vPudNET6U
# - [3]: Tangential Geometric Integral SPM, J(ij)
#           Link: https://www.youtube.com/watch?v=JRHnOsueic8
# - [4]: Streamline Geometric Integral SPM, Mx(ij) and My(ij)
#           Link: https://www.youtube.com/watch?v=BnPZjGCatcg
# - [5]: Solving the System of Equations
#           Link: https://www.youtube.com/watch?v=ep7vPzGYsbw
# - [6]: UIUC Airfoil Database: Download All Files using Python
#           Link: https://www.youtube.com/watch?v=nILo18DlqAo
# - [7]: Python code for downloading Selig airfoil DAT files
#           Link: http://www.joshtheengineer.com/2019/01/30/uiuc-airfoil-database-file-download/

import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path
from XFOIL import XFOIL
from COMPUTE_IJ_SPM import COMPUTE_IJ_SPM
from STREAMLINE_SPM import STREAMLINE_SPM
from COMPUTE_CIRCULATION import COMPUTE_CIRCULATION

# %% KNOWNS

# Airfoil loading flags
flagAirfoil = [0,                                                               # Create specified NACA airfoil in XFOIL
               1]                                                               # Load Selig-format airfoil from directory

# User-defined knowns
Vinf = 1                                                                        # Freestream velocity []
AoA  = 0                                                                        # Angle of attack [deg]
NACA = '2412'                                                                   # NACA airfoil to load [####]

# Convert angle of attack to radians
AoAR   = AoA*(np.pi/180)                                                        # Angle of attack [rad]

# Plotting flags
flagPlot = [1,      # Airfoil with panel normal vectors
            1,      # Geometry boundary pts, control pts, first panel, second panel
            1,      # Cp vectors at airfoil surface panels
            1,      # Pressure coefficient comparison (XFOIL vs. SPM)
            0,      # Airfoil streamlines
            0]      # Pressure coefficient contour

# %% XFOIL - CREATE/LOAD AIRFOIL

# PPAR menu options
PPAR = ['170',                                                                  # "Number of panel nodes"
        '4',                                                                    # "Panel bunching paramter"
        '1',                                                                    # "TE/LE panel density ratios"
        '1',                                                                    # "Refined area/LE panel density ratio"
        '1 1',                                                                  # "Top side refined area x/c limits"
        '1 1']                                                                  # "Bottom side refined area x/c limits"

# Call XFOIL function to obtain the following:
# - Airfoil coordinates
# - Pressure coefficient along airfoil surface
# - Lift, drag, and moment coefficients
xFoilResults = XFOIL(NACA, PPAR, AoA, flagAirfoil)                              # Get XFOIL results for prescribed airfoil

# Separate out results from XFOIL function results
afName  = xFoilResults[0]                                                       # Airfoil name
xFoilX  = xFoilResults[1]                                                       # X-coordinate for Cp result
xFoilY  = xFoilResults[2]                                                       # Y-coordinate for Cp result
xFoilCP = xFoilResults[3]                                                       # Pressure coefficient
XB      = xFoilResults[4]                                                       # Boundary point X-coordinate
YB      = xFoilResults[5]                                                       # Boundary point Y-coordinate
xFoilCL = xFoilResults[6]                                                       # Lift coefficient
xFoilCD = xFoilResults[7]                                                       # Drag coefficient
xFoilCM = xFoilResults[8]                                                       # Moment coefficient

# Number of boundary points and panels
numPts = len(XB)                                                                # Number of boundary points
numPan = numPts - 1                                                             # Number of panels (control points)

# %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

# Check for direction of points
edge = np.zeros(numPan)                                                         # Initialize edge value array
for i in range(numPan):                                                         # Loop over all panels
    edge[i] = (XB[i+1]-XB[i])*(YB[i+1]+YB[i])                                   # Compute edge values

sumEdge = np.sum(edge)                                                          # Sum all edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):                                                               # If panels are CCW
    XB = np.flipud(XB)                                                          # Flip the X-data array
    YB = np.flipud(YB)                                                          # Flip the Y-data array

# %% PANEL METHOD GEOMETRY - REF [1]
    
# Initialize variables
XC  = np.zeros(numPan)                                                          # Initialize control point X-coordinate
YC  = np.zeros(numPan)                                                          # Initialize control point Y-coordinate
S   = np.zeros(numPan)                                                          # Initialize panel length array
phi = np.zeros(numPan)                                                          # Initialize panel orientation angle array [deg]

# Find geometric quantities of the airfoil
for i in range(numPan):                                                         # Loop over all panels
    XC[i]   = 0.5*(XB[i]+XB[i+1])                                               # X-value of control point
    YC[i]   = 0.5*(YB[i]+YB[i+1])                                               # Y-value of control point
    dx      = XB[i+1]-XB[i]                                                     # Change in X between boundary points
    dy      = YB[i+1]-YB[i]                                                     # Change in Y between boundary points
    S[i]    = (dx**2 + dy**2)**0.5                                              # Length of the panel
    phi[i]  = math.atan2(dy,dx)                                                 # Angle of panel (positive X-axis to inside face)
    if (phi[i] < 0):                                                            # Make all panel angles positive [rad]
        phi[i] = phi[i] + 2*np.pi

# Compute angle of panel normal w.r.t. horizontal and include AoA
delta                = phi + (np.pi/2)                                          # Angle of panel normal [rad]
beta                 = delta - AoAR                                             # Angle of panel normal and AoA [rad]
beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi                           # Make all panel angles between 0 and 2pi [rad]

# %% COMPUTE SOURCE PANEL STRENGTHS - REF [5]

# Geometric integral (normal [I] and tangential [J])
# - Refs [2] and [3]
I, J = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S)                                        # Compute geometric integrals

# Populate A matrix
# - Simpler option: A = I + np.pi*np.eye(numPan,numPan)
A = np.zeros([numPan,numPan])                                                   # Initialize the A matrix
for i in range(numPan):                                                         # Loop over all i panels
    for j in range(numPan):                                                     # Loop over all j panels
        if (i == j):                                                            # If the panels are the same
            A[i,j] = np.pi                                                      # Set A equal to pi
        else:                                                                   # If panels are not the same
            A[i,j] = I[i,j]                                                     # Set A equal to geometric integral

# Populate b array
# - Simpler option: b = -Vinf*2*np.pi*np.cos(beta)
b = np.zeros(numPan)                                                            # Initialize the b array
for i in range(numPan):                                                         # Loop over all panels
    b[i] = -Vinf*2*np.pi*np.cos(beta[i])                                        # Compute RHS array

# Compute source panel strengths (lam) from system of equations
lam = np.linalg.solve(A,b)                                                      # Compute all source strength values

# Check the sum of the source strengths
# - This should be very close to zero for a closed polygon
sumLambda = sum(lam*S)                                                          # Check sum of source panel strengths
print("Sum of L: ",sumLambda)                                                   # Print sum of all source strengths

# %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

# Compute velocities
# - Simpler method: Vt = Vinf*np.sin(beta) + np.dot(J,lam)/(2*np.pi)
#                   Cp = 1 - (Vt/Vinf)**2
Vt = np.zeros(numPan)                                                           # Initialize tangential velocity array
Cp = np.zeros(numPan)                                                           # Initialize pressure coefficient array
for i in range(numPan):                                                         # Loop over all i panels
    addVal = 0                                                                  # Reset the summation value to zero
    for j in range(numPan):                                                     # Loop over all j panels
        addVal = addVal + (lam[j]/(2*np.pi))*J[i,j]                             # Sum all tangential source panel terms
    
    Vt[i] = Vinf*np.sin(beta[i]) + addVal                                       # Compute tangential velocity by adding uniform flow term
    Cp[i] = 1 - (Vt[i]/Vinf)**2                                                 # Compute pressure coefficient

# %% COMPUTE LIFT AND DRAG

# Compute normal and axial force coefficients
CN = -Cp*S*np.sin(beta)                                                         # Normal force coefficient []
CA = -Cp*S*np.cos(beta)                                                         # Axial force coefficient []

# Compute lift, drag, and moment coefficients
CL = sum(CN*np.cos(AoAR)) - sum(CA*np.sin(AoAR))                                # Decompose axial and normal to lift coefficient []
CD = sum(CN*np.sin(AoAR)) + sum(CA*np.cos(AoAR))                                # Decompose axial and normal to drag coefficient []
CM = sum(Cp*(XC-0.25)*S*np.cos(phi))                                            # Moment coefficient []

# Print the results to the Console
print("======= RESULTS =======")
print("Lift Coefficient (CL)")
print("\tSPM  : %2.8f" % CL)
print("\tXFOIL: %2.8f" % xFoilCL)
print("Drag Coefficient (CD)")
print("\tSPM  : %2.8f" % CD)
print("\tXFOIL: %2.8f" % xFoilCD)
print("Moment Coefficient (CM)")
print("\tSPM  : %2.8f" % CM)
print("\tXFOIL: %2.8f" % xFoilCM)

# %% COMPUTE STREAMLINES

if (flagPlot[4] == 1 or flagPlot[5] == 1):
    # Grid parameters
    nGridX   = 100                                                              # X-grid for streamlines and contours
    nGridY   = 100                                                              # Y-grid for streamlines and contours
    xVals    = [-0.5, 1.5]                                                      # X-grid extents [min, max]
    yVals    = [-0.5, 0.5]                                                      # Y-grid extents [min, max]
    
    # Streamline parameters
    slPct  = 40                                                                 # Percentage of streamlines of the grid
    Ysl    = np.linspace(yVals[0],yVals[1],int((slPct/100)*nGridY))             # Create array of Y streamline starting points
    Xsl    = xVals[0]*np.ones(len(Ysl))                                         # Create array of X streamline starting points
    XYsl   = np.vstack((Xsl.T,Ysl.T)).T                                         # Concatenate X and Y streamline starting points
    
    # Generate the grid points
    Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                              # X-values in evenly spaced grid
    Ygrid  = np.linspace(yVals[0],yVals[1],nGridY)                              # Y-values in evenly spaced grid
    XX, YY = np.meshgrid(Xgrid, Ygrid)                                          # Create meshgrid from X and Y grid arrays
    
    # Initialize velocities
    Vx     = np.zeros([nGridX,nGridY])                                          # Initialize X velocity matrix
    Vy     = np.zeros([nGridX,nGridY])                                          # Initialize Y velocity matrix
    
    # Path to figure out if grid point is inside polygon or not
    AF     = np.vstack((XB.T,YB.T)).T                                           # Concatenate XB and YB geometry points
    afPath = path.Path(AF)                                                      # Create a path for the geometry
    
    # Solve for grid point X and Y velocities
    for m in range(nGridX):                                                     # Loop over X grid points
        for n in range(nGridY):                                                 # Loop over Y grid points
            XP     = XX[m,n]                                                    # Current iteration's X grid point
            YP     = YY[m,n]                                                    # Current iteration's Y grid point
            Mx, My = STREAMLINE_SPM(XP,YP,XB,YB,phi,S)                          # Compute Mx and My geometric integrals
            
            # Check if grid points are in object
            # - If they are, assign a velocity of zero
            if afPath.contains_points([(XP,YP)]):                               # If (XP,YP) is in the body
                Vx[m,n] = 0                                                     # X-velocity is zero
                Vy[m,n] = 0                                                     # Y-velocity is zero
            else:
                Vx[m,n] = Vinf*np.cos(AoAR) + sum(lam*Mx/(2*np.pi))             # Set X-velocity
                Vy[m,n] = Vinf*np.sin(AoAR) + sum(lam*My/(2*np.pi))             # Set Y-velocity
    
    # Compute grid point velocity magnitude and pressure coefficient
    Vxy  = np.sqrt(Vx**2 + Vy**2)                                               # Compute magnitude of velocity vector
    CpXY = 1 - (Vxy/Vinf)**2                                                    # Pressure coefficient []

# %% CIRCULATION AND SOURCE STRENGTH CHECK

if (flagPlot[4] == 1 or flagPlot[5] == 1):
    # Compute circulation
    a    = 0.75                                                                 # Ellipse horizontal half-length
    b    = 0.25                                                                 # Ellipse vertical half-length
    x0   = 0.5                                                                  # Ellipse center X-coordinate
    y0   = 0                                                                    # Ellipse center Y-coordinate
    numT = 5000                                                                 # Number of points on ellipse
    Circulation, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(a,b,x0,y0,              # Computer circulation around ellipse
                                                        numT,Vx,Vy,Xgrid,Ygrid)
    
    # Print values to Console
    print("Sum of L   : %2.8f" % sumLambda)                                     # Print sum of source strengths
    print("Circulation: %2.8f" % Circulation)                                   # Print circulation

# %% PLOTTING

# FIGURE: Airfoil with panel normal vectors
if (flagPlot[0] == 1):
    fig = plt.figure(1)                                                         # Create the figure
    plt.cla()                                                                   # Clear the axes
    plt.fill(XB,YB,'k')                                                         # Plot the airfoil
    X = np.zeros(2)                                                             # Initialize 'X'
    Y = np.zeros(2)                                                             # Initialize 'Y'
    for i in range(numPan):                                                     # Loop over all panels
        X[0] = XC[i]                                                            # Set X start of panel orientation vector
        X[1] = XC[i] + S[i]*np.cos(delta[i])                                    # Set X end of panel orientation vector
        Y[0] = YC[i]                                                            # Set Y start of panel orientation vector
        Y[1] = YC[i] + S[i]*np.sin(delta[i])                                    # Set Y end of panel orientation vector
        if (i == 0):                                                            # If it's the first panel index
            plt.plot(X,Y,'b-',label='First Panel')                              # Plot normal vector for first panel
        elif (i == 1):                                                          # If it's the second panel index
            plt.plot(X,Y,'g-',label='Second Panel')                             # Plot normal vector for second panel
        else:                                                                   # If it's neither the first nor second panel index
            plt.plot(X,Y,'r-')                                                  # Plot normal vector for all other panels
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Panel Geometry')                                                 # Set title
    plt.axis('equal')                                                           # Set axes equal
    plt.legend()                                                                # Display legend
    plt.show()                                                                  # Display plot

# FIGURE: Geometry with the following indicated:
# - Boundary points, control points, first panel, second panel
if (flagPlot[1] == 1):
    fig = plt.figure(2)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.plot(XB,YB,'k-')                                                        # Plot airfoil panels
    plt.plot([XB[0], XB[1]],[YB[0], YB[1]],'b-',label='First Panel')            # Plot first panel
    plt.plot([XB[1], XB[2]],[YB[1], YB[2]],'g-',label='Second Panel')           # Plot second panel
    plt.plot(XB,YB,'ko',markerfacecolor='k',label='Boundary Pts')               # Plot boundary points (black circles)
    plt.plot(XC,YC,'ko',markerfacecolor='r',label='Control Pts')                # Plot control points (red circles)
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.axis('equal')                                                           # Set axes equal
    plt.legend()                                                                # Display legend
    plt.show()                                                                  # Display plot
    
# FIGURE: Cp vectors at airfoil control points
if (flagPlot[2] == 1):
    fig = plt.figure(3)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    Cps = np.absolute(Cp*0.15)                                                  # Scale and make positive all Cp values
    X = np.zeros(2)                                                             # Initialize X values
    Y = np.zeros(2)                                                             # Initialize Y values
    for i in range(len(Cps)):                                                   # Loop over all panels
        X[0] = XC[i]                                                            # Control point X-coordinate
        X[1] = XC[i] + Cps[i]*np.cos(delta[i])                                  # Ending X-value based on Cp magnitude
        Y[0] = YC[i]                                                            # Control point Y-coordinate
        Y[1] = YC[i] + Cps[i]*np.sin(delta[i])                                  # Ending Y-value based on Cp magnitude
        
        if (Cp[i] < 0):                                                         # If pressure coefficient is negative
            plt.plot(X,Y,'r-')                                                  # Plot as a red line
        elif (Cp[i] >= 0):                                                      # If pressure coefficient is zero or positive
            plt.plot(X,Y,'b-')                                                  # Plot as a blue line
    plt.fill(XB,YB,'k')                                                         # Plot the airfoil as black polygon
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.gca().set_aspect('equal')                                               # Set aspect ratio equal
    plt.show()                                                                  # Show the plot

# FIGURE: Cp vectors at airfoil surface panels
if (flagPlot[3] == 1):
    fig = plt.figure(4)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    midIndX = int(np.floor(len(xFoilCP)/2))                                     # Airfoil middle index for XFOIL data
    midIndS = int(np.floor(len(Cp)/2))                                          # Airfoil middle index for SPM data
    plt.plot(xFoilX[0:midIndX],xFoilCP[0:midIndX],                              # Plot Cp for upper surface of airfoil from XFOIL
             'b-',label='Upper')
    plt.plot(xFoilX[midIndX+1:len(xFoilX)],xFoilCP[midIndX+1:len(xFoilX)],      # Plot Cp for lower surface of airfoil from XFOIL
             'r-',label='Lower')
    plt.plot(XC[0:midIndS],Cp[0:midIndS],                                       # Plot Cp for upper surface of airfoil from panel method
             'ks',markerfacecolor='r')
    plt.plot(XC[midIndS+1:len(XC)],Cp[midIndS+1:len(XC)],                       # Plot Cp for lower surface of airfoil from panel method
             'ks',markerfacecolor='b')
    plt.xlim(0,1)                                                               # Set X-limits
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Pressure Coefficient')                                           # Set title
    plt.show()                                                                  # Display plot
    plt.legend()                                                                # Display legend
    plt.gca().invert_yaxis()                                                    # Invert Cp (Y) axis

# FIGURE: Airfoil streamlines
if (flagPlot[4] == 1):
    fig = plt.figure(5)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.streamplot(XX,YY,Vx,Vy, linewidth=0.5, density=10, color='r',           # Plot streamlines
                   arrowstyle='-', start_points=XYsl)
    plt.fill(XB,YB,'k')                                                         # Plot airfoil as black polygon
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()                                                                  # Display plot

# FIGURE: Pressure coefficient contour
if (flagPlot[5] == 1):
    fig = plt.figure(6)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.contourf(XX,YY,CpXY,500,cmap='jet')                                     # Plot contour
    plt.fill(XB,YB,'k')                                                         # Plot airfoil as black polygon
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()                                                                  # Display plot

