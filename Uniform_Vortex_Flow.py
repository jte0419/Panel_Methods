# ELEMENTARY FLOW - UNIFORM FLOW + VORTEX FLOW
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started: 02/19/19
# Updated: 02/19/19 - Transferred from MATLAB to Python
#                   - Works as expected

import numpy as np
import matplotlib.pyplot as plt
from COMPUTE_CIRCULATION import COMPUTE_CIRCULATION

# %% KNOWNS

Vinf  = 1                                                                       # Freestream velocity [arb]
alpha = 0                                                                       # Angle of attack [deg]
gamma = 10                                                                      # Vortex strength (+: CW, -: CCW)
X0    = 0                                                                       # Vortex X coordinate
Y0    = 0                                                                       # Vortex Y coordinate

# %% CALCULATIONS

# Create grid
numX   = 50                                                                     # Number of X points
numY   = 50                                                                     # Number of Y points
X      = np.linspace(-10,10,numX)                                               # X-point array
Y      = np.linspace(-10,10,numY)                                               # Y-point array
XX, YY = np.meshgrid(X,Y)                                                       # Create the meshgrid

# Solve for velocities
Vx = np.zeros([numX,numY])                                                      # Initialize
Vy = np.zeros([numX,numY])                                                      # Initialize
r  = np.zeros([numX,numY])                                                      # Initialize
for i in range(numX):                                                           # Loop over X points
    for j in range(numY):                                                       # Loop over Y points
        x       = XX[i,j]                                                       # X coordinate
        y       = YY[i,j]                                                       # Y coordinate
        dx      = x - X0                                                        # X distance from vortex
        dy      = y - Y0                                                        # Y distance from vortex
        r       = np.sqrt(dx**2 + dy**2)                                        # Distance from vortex
        Vx[i,j] = Vinf*np.cos(alpha*(np.pi/180)) +  (gamma*dy)/(2*np.pi*r**2)   # Compute X velocity component
        Vy[i,j] = Vinf*np.sin(alpha*(np.pi/180)) +  (-gamma*dx)/(2*np.pi*r**2)  # Compute Y velocity component

# %% COMPUTE CIRCULATIONS

# First circulation ellipse
a    = 5                                                                        # Horizontal axis half-length
b    = 5                                                                        # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 0                                                                        # Ellipse center Y coordinate
numT = 5000                                                                     # Number of points along ellipse
Gamma1, xC1, yC1, VxC1, VyC1 = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)    # Call circulation calculation 
print("Circulation 1: ", Gamma1)                                                # Display circulation 1 result

# Second circulation ellipse
a    = 1.5                                                                      # Horizontal axis half-length
b    = 1.5                                                                      # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 3                                                                        # Ellipse center Y coordinate
numT = 5000                                                                     # Number of points along ellipse
Gamma2, xC2, yC2, VxC2, VyC2 = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)    # Call circulation calculation
print("Circulation 2: ", Gamma2)                                                # Display circulation 2 result

# %% PLOTTING

# Streamline starting points
numSL = 20                                                                      # Number of streamlines
Xsl  = -10*np.ones(numSL)                                                       # Streamline starting X coordinates
Ysl  = np.linspace(-10,10,numSL)                                                # Streamline starting Y coordinates
XYsl = np.vstack((Xsl.T,Ysl.T)).T                                               # Concatenate X and Y streamline starting points

# Plot quiver and streamlines
fig = plt.figure(1)                                                             # Create the figure
plt.cla()                                                                       # Get ready for plotting
plt.quiver(X,Y,Vx,Vy)                                                           # Plot velocity vectors
plt.streamplot(XX,YY,Vx,Vy, linewidth=0.5, density=10, color='r', arrowstyle='-', start_points=XYsl)    # Plot streamlines
line_1, = plt.plot(xC1,yC1,'b-')                                                # Plot first ellipse
line_2, = plt.plot(xC2,yC2,'m-')                                                # Plot second ellipse
plt.title('Uniform + Vortex Flow')                                              # Set title
plt.xlabel('X-Axis')                                                            # Set X-label
plt.ylabel('Y-Axis')                                                            # Set Y-label
plt.gca().set_aspect('equal')                                                   # Set axes equal
plt.xlim([-6, 6])                                                               # Set X-limits
plt.ylim([-6, 6])                                                               # Set Y-limits
leg1Str = '$\Gamma$ = ' + "{:.3f}".format(Gamma1)
leg2Str = '$\Gamma$ = ' + "{:.3f}".format(Gamma2)
plt.legend([line_1, line_2],[leg1Str,leg2Str], loc='lower right', facecolor='white')
plt.show()                                                                      # Display plot

## Save the figure
#saveFlnm = 'Uniform_Vortex_30'
#savePth = 'C:/Users/Josh/Documents/Latex/YouTube/Panel_Methods/Figures_Python/' + saveFlnm + '.pdf'
#plt.savefig(savePth,bbox_inches='tight',facecolor='w',edgecolor='w')
