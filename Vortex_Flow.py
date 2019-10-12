# ELEMENTARY FLOW - VORTEX FLOW
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

gamma = 2                                                                       # Vortex strength (+: CW, -: CCW)
X0    = 0                                                                       # Vortex X coordinate
Y0    = 0                                                                       # Vortex Y coordinate

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
V  = np.zeros([numX,numY])                                                      # Initialize velocity
Vt = np.zeros([numX,numY])                                                      # Initialize tangential velocity component
Vr = np.zeros([numX,numY])                                                      # Initialize radial velocity component
r  = np.zeros([numX,numY])                                                      # Initialize radius
for i in range(numX):                                                           # Loop over X points
    for j in range(numY):                                                       # Loop over Y points
        x       = XX[i,j]                                                       # X coordinate
        y       = YY[i,j]                                                       # Y coordinate
        dx      = x - X0                                                        # X distance from vortex
        dy      = y - Y0                                                        # Y distance from vortex
        r       = np.sqrt(dx**2 + dy**2)                                        # Distance from vortex
        Vx[i,j] = (gamma*dy)/(2*np.pi*r**2)                                     # Compute X velocity component
        Vy[i,j] = (-gamma*dx)/(2*np.pi*r**2)                                    # Compute Y velocity component
        V[i,j]  = np.sqrt(Vx[i,j]**2 + Vy[i,j]**2)                              # Compute velocity
        Vt[i,j] = -gamma/(2*np.pi*r)                                            # Compute tangential velocity component
        Vr[i,j] = 0                                                             # Compute radial velocity component

# %% COMPUTE CIRCULATIONS

a    = 2                                                                        # Horizontal axis half-length
b    = 2                                                                        # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 0                                                                        # Ellipse center Y coordinate
numT = 50                                                                       # Number of points along ellipse
Gamma, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)         # Call circulation calculation
print("Circulation: ", Gamma)                                                   # Display circulation result
        
#a    = 0.5                                                                      # Horizontal axis half-length
#b    = 0.5                                                                      # Vertical axis half-length
#x0   = 0                                                                        # Ellipse center X coordinate
#y0   = 2                                                                        # Ellipse center Y coordinate
#numT = 50                                                                       # Number of points along ellipse
#Gamma, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)         # Call circulation calculation
#print("Circulation: ", Gamma)                                                   # Display circulation result

# %% PLOTTING

# Plot velocity vectors
fig = plt.figure(1)                                                             # Create figure
plt.cla()                                                                       # Get ready for plotting
plt.quiver(X,Y,Vx,Vy)                                                           # Plot velocity vectors
plt.plot(xC,yC,'b-')                                                            # Plot ellipse
plt.title('Vortex Flow')                                                        # Set title
plt.xlabel('X-Axis')                                                            # Set X-label
plt.ylabel('Y-Axis')                                                            # Set Y-label
plt.xlim([-3, 3])                                                               # Set X-limits
plt.ylim([-3, 3])                                                               # Set Y-limits
plt.gca().set_aspect('equal')                                                   # Set axes equal
plt.show()                                                                      # Display plot

## Save the figure
#saveFlnm = 'Vortex_Flow_Pos'
#savePth = 'C:/Users/Josh/Documents/Latex/YouTube/Panel_Methods/Figures_Python/' + saveFlnm + '.pdf'
#plt.savefig(savePth,bbox_inches='tight',facecolor='w',edgecolor='w')
