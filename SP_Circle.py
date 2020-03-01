# SOURCE PANEL METHOD - CIRCLE
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 01/01/19
# Updated   : 01/01/19 - Started code in MATLAB
#                      - Works as expected
#             02/03/19 - Transferred from MATLAB to Python
#                      - Works as expected
# Notes     : This code is not optimized, but is instead written in such a way
#             that it is easy to follow along with my YouTube video derivations
# 
# Functions Needed:
# - COMPUTE_IJ_SPM.py
# - STREAMLINE_SPM.py
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

import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path
from COMPUTE_IJ_SPM import COMPUTE_IJ_SPM
from STREAMLINE_SPM import STREAMLINE_SPM

# %% KNOWNS

# User-defined knowns
Vinf = 1                                                                        # Freestream velocity
AoA  = 0                                                                        # Angle of attack [deg]
numB = 9                                                                        # Number of boundary points (including endpoint)
tO   = (360/(numB-1))/2                                                         # Boundary point angle offset [deg]
AoAR = AoA*(np.pi/180)                                                          # Convert AoA to radians [rad]

# Plotting flags
flagPlot = [1,      # Shape polygon with panel normal vectors
            1,      # Geometry boundary pts, control pts, first panel, second panel
            1,      # Analytical and SPM pressure coefficient plot
            1,      # Streamline plot
            1]      # Pressure coefficient contour plot

# %% CREATE CIRCLE BOUNDARY POINTS

# Angles used to compute boundary points
theta = np.linspace(0,360,numB)                                                 # Create angles for computing boundary point locations [deg]
theta = theta + tO                                                              # Add panel angle offset [deg]
theta = theta*(np.pi/180)                                                       # Convert from degrees to radians [rad]

# Boundary points
XB = np.cos(theta)                                                              # Compute boundary point X-coordinate [radius of 1]
YB = np.sin(theta)                                                              # Compute boundary point Y-coordinate [radius of 1]

# Number of panels
numPan = len(XB)-1                                                              # Number of panels (control points)

# %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

# Check for direction of points
edge = np.zeros(numPan)                                                         # Initialize edge value array
for i in range(numPan):                                                         # Loop over all panels
    edge[i] = (XB[i+1]-XB[i])*(YB[i+1]+YB[i])                                   # Compute edge values

sumEdge = np.sum(edge)                                                          # Sum of all edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):                                                               # If panels are CCW
    XB = np.flipud(XB)                                                          # Flip the X-data array
    YB = np.flipud(YB)                                                          # Flip the Y-data array

# %% PANEL METHOD GEOMETRY - REF [1]

# Initialize variables
XC  = np.zeros(numPan)                                                          # Initialize control point X-coordinate
YC  = np.zeros(numPan)                                                          # Initialize control point Y-coordinate
S   = np.zeros(numPan)                                                          # Intialize panel length array
phi = np.zeros(numPan)                                                          # Initialize panel orientation angle array

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
print("Sum of L: ",sum(lam*S))                                                  # Print sum of all source strengths

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

# Analytical angles and pressure coefficients
analyticTheta = np.linspace(0,2*np.pi,200)                                      # Analytical theta angles [rad]
analyticCP    = 1 - 4*np.sin(analyticTheta)**2                                  # Analytical pressure coefficient []

# %% COMPUTE LIFT AND DRAG

# Compute normal and axial force coefficients
CN = -Cp*S*np.sin(beta)                                                         # Normal force coefficient []
CA = -Cp*S*np.cos(beta)                                                         # Axial force coefficient []

# Compute lift and drag coefficients
CL = sum(CN*np.cos(AoAR)) - sum(CA*np.sin(AoAR))                                # Decompose axial and normal to lift coefficient []
CD = sum(CN*np.sin(AoAR)) + sum(CA*np.cos(AoAR))                                # Decompose axial and normal to drag coefficient []

print("CL      : ",CL)                                                          # Display lift coefficient (should be zero)
print("CD      : ",CD)                                                          # Display drag coefficient (should be zero)

# %% COMPUTE STREAMLINES - REF [4]

if (flagPlot[3] == 1 or flagPlot[4] == 1):
    # Grid parameters
    nGridX   = 100                                                              # X-grid for streamlines and contours
    nGridY   = 100                                                              # Y-grid for streamlines and contours
    xVals    = [-1.5, 1.5]                                                      # X-grid extents [min, max]
    yVals    = [-1.5, 1.5]                                                      # Y-grid extents [min, max]
    
    # Streamline parameters
    slPct  = 30                                                                 # Percentage of streamlines of the grid
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
    for m in range(nGridX):                                                     # Loop over X-grid points
        for n in range(nGridY):                                                 # Loop over Y-grid points
            XP     = XX[m,n]                                                    # Isolate X point
            YP     = YY[m,n]                                                    # Isolate Y point
            Mx, My = STREAMLINE_SPM(XP,YP,XB,YB,phi,S)                          # Compute streamline Mx and My values (Ref [4])
            
            # Check if grid points are in object
            # - If they are, assign a velocity of zero
            if afPath.contains_points([(XP,YP)]):                               # If (XP,YP) is in the polygon body
                Vx[m,n] = 0                                                     # X-velocity is zero
                Vy[m,n] = 0                                                     # Y-velocity is zero
            else:                                                               # If (XP,YP) is not in the polygon body
                Vx[m,n] = Vinf*np.cos(AoAR) + sum(lam*Mx/(2*np.pi))             # Compute X-velocity
                Vy[m,n] = Vinf*np.sin(AoAR) + sum(lam*My/(2*np.pi))             # Compute Y-velocity
            
    # Compute grid point velocity magnitude and pressure coefficient
    Vxy  = np.sqrt(Vx**2 + Vy**2)                                               # Compute magnitude of velocity vector
    CpXY = 1 - (Vxy/Vinf)**2                                                    # Pressure coefficient []

# %% PLOTTING

# FIGURE: Shape polygon with panel normal vectors
if (flagPlot[0] == 1):
    angCirc = np.linspace(0,2*np.pi,1000)                                       # Angles for "perfect" circle
    xCirc = np.cos(angCirc)                                                     # "Perfect" circle X values
    yCirc = np.sin(angCirc)                                                     # "Perfect" circle Y values
    fig = plt.figure(1)                                                         # Create figure
    plt.cla()                                                                   # Clear the axes
    plt.plot(xCirc,yCirc,'k--')                                                 # Plot the circle that polygon is approximating
    plt.fill(XB,YB,'k')                                                         # Plot the paneled circle
    X = np.zeros(2)                                                             # Initialize 'X'
    Y = np.zeros(2)                                                             # Initialize 'Y'
    for i in range(numPan):                                                     # Loop over all panels
        X[0] = XC[i]                                                            # Set X start of panel orientation vector
        X[1] = XC[i] + S[i]*np.cos(delta[i])                                    # Set X end of panel orientation vector
        Y[0] = YC[i]                                                            # Set Y start of panel orientation vector
        Y[1] = YC[i] + S[i]*np.sin(delta[i])                                    # Set Y end of panel orientation vector
        if (i == 0):                                                            # If it's the first panel index
            plt.plot(X,Y,'b-',label='First Panel')                              # Plot the first panel
        elif (i == 1):                                                          # If it's the second panel index
            plt.plot(X,Y,'g-',label='Second Panel')                             # Plot the second panel
        else:                                                                   # If it's neither the first nor second panel index
            plt.plot(X,Y,'r-')                                                  # Plot the rest of the panels
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Panel Geometry')                                                 # Set title
    plt.axis('equal')                                                           # Set axes equal
    plt.legend()                                                                # Show legend
    plt.show()                                                                  # Display plot

# FIGURE: Geometry with the following indicated:
# - Boundary points, control points, first panel, second panel
if (flagPlot[1] == 1):
    fig = plt.figure(2)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.plot(XB,YB,'k-',label='Panels')                                         # Plot polygon
    plt.plot([XB[0], XB[1]],[YB[0], YB[1]],'b-',label='First Panel')            # Plot first panel
    plt.plot([XB[1], XB[2]],[YB[1], YB[2]],'g-',label='Second Panel')           # Plot second panel
    plt.plot(XB,YB,'ko',markerfacecolor='k',label='Boundary Points')            # Plot boundary points
    plt.plot(XC,YC,'ko',markerfacecolor='r',label='Control Points')             # Plot control points
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Panel Geometry 2')                                               # Set title
    plt.axis('equal')                                                           # Set axes equal
    plt.legend()                                                                # Show legend
    plt.show()                                                                  # Display plot

# FIGURE: Analytical and SPM pressure coefficient
if (flagPlot[2] == 1):
    fig = plt.figure(3)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.plot(analyticTheta*(180/np.pi),analyticCP,'b-',label='Analytical')      # Plot analytical pressure coefficient
    plt.plot(beta*(180/np.pi),Cp,'ks',markerfacecolor='r',label='SPM')          # Plot panel method pressure coefficient
    plt.xlabel('Angle [deg]')                                                   # Set X-label
    plt.ylabel('Pressure Coefficient')                                          # Set Y-label
    plt.title('Pressure Coefficient Comparison')                                # Set title
    plt.xlim(0, 360)                                                            # Set X-limits
    plt.ylim(-3.5, 1.5)                                                         # Set Y-limits
    plt.legend()                                                                # Show legend
    plt.show()                                                                  # Display plot

# FIGURE: Streamlines
if (flagPlot[3] == 1):
    fig = plt.figure(4)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.streamplot(XX,YY,Vx,Vy, linewidth=0.5, density=10, color='r', 
                   arrowstyle='-', start_points=XYsl)                           # Plot the streamlines
    plt.fill(XB,YB,'k')                                                         # Plot polygon
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Streamline')                                                     # Set title
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()                                                                  # Display plot

# FIGURE: Pressure coefficient contours
if (flagPlot[4] == 1):
    fig = plt.figure(5)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.contourf(XX,YY,CpXY,500,cmap='jet')                                     # Plot contour
    plt.fill(XB,YB,'k')                                                         # Plot polygon
    plt.xlabel('X-Axis')                                                        # Set X-label
    plt.ylabel('Y-Axis')                                                        # Set Y-label
    plt.title('Pressure Coefficient Contours')                                  # Set title
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()                                                                  # Display plot
