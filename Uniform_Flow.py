# ELEMENTARY FLOW - UNIFORM FLOW
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

# %% CALCULATIONS

# Create grid
numX   = 10                                                                     # Number of X points
numY   = 10                                                                     # Number of Y points
X      = np.linspace(-10,10,numX)                                               # X-point array
Y      = np.linspace(-10,10,numY)                                               # Y-point array
XX, YY = np.meshgrid(X,Y)                                                       # Create the meshgrid

# Solve for velocities
Vx = np.zeros([numX,numY])                                                      # Initialize X velocity component
Vy = np.zeros([numX,numY])                                                      # Initialize Y velocity component
for i in range(numX):                                                           # Loop over X points
    for j in range(numY):                                                       # Loop over Y points
        Vx[i,j] = Vinf*np.cos(alpha*(np.pi/180))                                # Compute X velocity components
        Vy[i,j] = Vinf*np.sin(alpha*(np.pi/180))                                # Compute Y velocity components

# %% COMPUTE CIRCULATION

a    = 5                                                                        # Horizontal axis half-length
b    = 5                                                                        # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 0                                                                        # Ellipse center Y coordinate
numT = 100                                                                      # Number of points along ellipse
Gamma, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)         # Call circulation calculation
print("Circulation: ", Gamma)                                                   # Display circulation result

# %% PLOTTING

# Plot quiver
fig = plt.figure(1)                                                             # Create figure
plt.cla()                                                                       # Get ready for plotting
plt.quiver(X,Y,Vx,Vy)                                                           # Plot velocity vectors
plt.plot(xC,yC,'b-')                                                            # Plot ellipse
plt.title('Uniform Flow')                                                       # Set title
plt.xlabel('X-Axis')                                                            # Set X-label
plt.ylabel('Y-Axis')                                                            # Set Y-label
plt.gca().set_aspect('equal')                                                   # Set axes equal
plt.show()                                                                      # Display plot

## Save the figure
#saveFlnm = 'Uniform_Flow_30deg'
#savePth = 'C:/Users/Josh/Documents/Latex/YouTube/Panel_Methods/Figures_Python/' + saveFlnm + '.pdf'
#plt.savefig(savePth,bbox_inches='tight',facecolor='w',edgecolor='w')
