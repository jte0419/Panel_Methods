# ELEMENTARY FLOW - SOURCE/SINK FLOW
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

lmbda = 1                                                                       # Source/sink strength (+: Source, -: Sink)
X0    = 0                                                                       # Source/sink X coordinate
Y0    = 0                                                                       # Source/sink Y coordinate

# %% CALCULATIONS

# Create grid
numX   = 100                                                                    # Number of X points
numY   = 100                                                                    # Number of Y points
X      = np.linspace(-10,10,numX)                                               # X-point array
Y      = np.linspace(-10,10,numY)                                               # Y-point array
XX, YY = np.meshgrid(X,Y)                                                       # Create the meshgrid

# Solve for velocities
Vx = np.zeros([numX,numY])                                                      # Initialize X velocity component
Vy = np.zeros([numX,numY])                                                      # Initialize Y velocity component
V  = np.zeros([numX,numY])                                                      # Initialize velocity magnitude
Vr = np.zeros([numX,numY])                                                      # Initialize radial velocity component
r  = np.zeros([numX,numY])                                                      # Initialize radius
for i in range(numX):                                                           # Loop over X points
    for j in range(numY):                                                       # Loop over Y points
        x       = XX[i,j]                                                       # X coordinate
        y       = YY[i,j]                                                       # Y coordinate
        dx      = x - X0                                                        # X distance from source/sink
        dy      = y - Y0                                                        # Y distance from source/sink
        r       = np.sqrt(dx**2 + dy**2)                                        # Distance from source/sink
        Vx[i,j] = (lmbda*dx)/(2*np.pi*r**2)                                     # Compute X velocity component
        Vy[i,j] = (lmbda*dy)/(2*np.pi*r**2)                                     # Compute Y velocity component
        V[i,j]  = np.sqrt(Vx[i,j]**2 + Vy[i,j]**2)                              # Compute velocity
        Vr[i,j] = lmbda/(2*np.pi*r)                                             # Compute radial velocity component

# %% COMPUTE CIRCULATIONS

a    = 1.5                                                                      # Horizontal axis half-length
b    = 1.5                                                                      # Vertical axis half-length
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
plt.title('Source/Sink Flow')                                                   # Set title
plt.xlabel('X-Axis')                                                            # Set X-label
plt.ylabel('Y-Axis')                                                            # Set Y-label
plt.xlim([-2, 2])                                                               # Set X-limits
plt.ylim([-2, 2])                                                               # Set Y-limits
plt.gca().set_aspect('equal')                                                   # Set axes equal
plt.show()                                                                      # Display plot

## Save the figure
#saveFlnm = 'Sink_Flow'
#savePth = 'C:/Users/Josh/Documents/Latex/YouTube/Panel_Methods/Figures_Python/' + saveFlnm + '.pdf'
#plt.savefig(savePth,bbox_inches='tight',facecolor='w',edgecolor='w')
