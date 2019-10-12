# CIRCULATION
# Written by: JoshTheEngineer
# Started: 03/03/19
# Updated: 03/03/19 - Started code
#                   - Works as expected

import numpy as np
import matplotlib.pyplot as plt
from COMPUTE_CIRCULATION import COMPUTE_CIRCULATION

# %% KNOWNS

N = 5;                                                                          # Number of X and Y grid points

X  = np.linspace(0,N-1,N);                                                      # X grid points
Y  = np.linspace(0,N-1,N);                                                      # Y grid points
Vx = np.random.rand(N,N);                                                       # Random X velocities
Vy = np.random.rand(N,N);                                                       # Random Y velocities

XX, YY = np.meshgrid(X,Y)                                                       # Create the meshgrid

# Create ellipse
a    = 1.75;                                                                    # Horizontal half-length
b    = 1.25;                                                                    # Vertical half-length
x0   = 2;                                                                       # Center X-point
y0   = 2;                                                                       # Center Y-point
numT = 50;                                                                      # Number of ellipse points

Gamma, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y)         # Compute circulation
print("Circulation: ", Gamma)                                                   # Print circulation

# %% PLOTTING

# Plot velocity vectors
fig = plt.figure(1)                                                             # Create figure
plt.cla()                                                                       # Get ready for plotting
plt.quiver(XX,YY,Vx,Vy,color='r')                                               # Velocity vectors
plt.plot(XX,YY,'k.')                                                            # Grid points
plt.quiver(xC,yC,VxC,VyC,color='b')                                             # Ellipse velocity vectors
plt.plot(xC,yC,'b.');                                                           # Ellipse points
plt.title('Circulation')                                                        # Title
plt.xlabel('X-Axis')                                                            # X-axis label
plt.ylabel('Y-Axis')                                                            # Y-axis label
plt.xlim([-1, 5])                                                               # X-axis limits
plt.ylim([-1, 5])                                                               # Y-axis limits
plt.gca().set_aspect('equal')                                                   # Set axes equal
plt.show()                                                                      # Display plot
