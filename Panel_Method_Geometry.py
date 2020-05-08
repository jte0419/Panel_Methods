# PANEL METHOD GEOMETRY
# Written by: JoshTheEngineer
# Started: 01/24/19
# Updated: 01/24/19 - Started code in MATLAB
#                   - Works as expected
#          02/03/19 - Transferred code to Python
#                   - Works as expected

import numpy as np
import math as m
import matplotlib.pyplot as plt
from LOAD_AIRFOIL_SELIG import LOAD_AIRFOIL_SELIG

# %% CREATE/LOAD GEOMETRY

# Load circle or airfoil
#load = 'Airfoil'                                                               # Use the airfoil geometry
load = 'Circle'                                                                 # Use the circle geometry
AoA  = 0                                                                        # Angle of attack [deg]

if (load == 'Airfoil'):                                                         # If airfoil is to be loaded
    fileName = 'goe623'                                                         # Specify the airfoil filename
    XB, YB = LOAD_AIRFOIL_SELIG(fileName)                                       # Load X and Y data points of airfoil
elif (load == 'Circle'):                                                        # If circle is to be used
    numB = 9                                                                    # Number of boundary points
    tO   = 22.5                                                                 # Angle offset [deg]
    
    # Angles used to compute boundary points
    theta = np.linspace(0,360,numB)                                             # Angles to compute boundary points [deg]
    theta = theta + tO                                                          # Add angle offset [deg]
    theta = theta*(np.pi/180)                                                   # Convert angle to radians [rad]

    # Boundary points
    XB = np.cos(theta)                                                          # X value of boundary points
    YB = np.sin(theta)                                                          # Y value of boundary points

# Number of panels
numPan =len(XB)-1                                                               # Number of panels

# Check for direction of points
edge = np.zeros(numPan)                                                         # Initialize edge check value
for i in range(numPan):                                                         # Loop over all panels
    edge[i] = (XB[i+1]-XB[i])*(YB[i+1]+YB[i])                                   # Compute edge value for each panel

sumEdge = np.sum(edge)                                                          # Sum all panel edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):                                                               # If sum is negative
    print('Points are counter-clockwise.  Flipping.\n')                         # Display message in console
    XB = np.flipud(XB)                                                          # Flip the X boundary points array
    YB = np.flipud(YB)                                                          # Flip the Y boundary points array
elif (sumEdge > 0):                                                             # If sum is positive
    print('Points are clockwise.  Not flipping.\n')                             # Do nothing, display message in consolve
    
# %% COMPUTE GEOMETRIC VARIABLES

# Initialize variables
XC  = np.zeros(numPan)                                                          # Initialize X control points
YC  = np.zeros(numPan)                                                          # Initialize Y control points
S   = np.zeros(numPan)                                                          # Initialize panel lengths
phi = np.zeros(numPan)                                                          # Initialize panel orientation angles

# Find geometric quantities of the airfoil
for i in range(numPan):                                                         # Loop over all panels
    XC[i]   = 0.5*(XB[i]+XB[i+1])                                               # X control point coordinate
    YC[i]   = 0.5*(YB[i]+YB[i+1])                                               # Y control point coordinate
    dx      = XB[i+1]-XB[i]                                                     # Panel X length
    dy      = YB[i+1]-YB[i]                                                     # Panel Y length
    S[i]    = (dx**2 + dy**2)**0.5                                              # Panel length
    phi[i]  = m.atan2(dy,dx)                                                    # Panel orientation angle [rad]
    if (phi[i] < 0):                                                            # If panel orientation is negative
        phi[i] = phi[i] + 2*np.pi                                               # Add 2pi to the panel angle

# Compute angle of panel normal w.r.t. horizontal and include AoA
delta = phi + (np.pi/2)                                                         # Compute panel normal angle [rad]
beta  = delta - (AoA*(np.pi/180))                                               # Angle between freestream and panel normal [rad]

# %% PLOTTING

# Dashed circle defined
T = np.linspace(0,2*np.pi,1000)                                                 # Angle array to compute circle
x = np.cos(T)                                                                   # Circle X points
y = np.sin(T)                                                                   # Circle Y points

# Plot the paneled geometry
fig = plt.figure(1)                                                             # Create figure
plt.cla()                                                                       # Get ready for plotting
if (load == 'Circle'):                                                          # If circle is selected
    plt.plot(x,y,'k--')                                                         # Plot actual circle outline
plt.fill(XB,YB,'k')                                                             # Plot polygon (circle or airfoil)
X = np.zeros(2)                                                                 # Initialize panel X variable
Y = np.zeros(2)                                                                 # Initialize panel Y variable
for i in range(numPan):                                                         # Loop over all panels
    X[0] = XC[i]                                                                # Panel starting X point
    X[1] = XC[i] + S[i]*np.cos(delta[i])                                        # Panel ending X point
    Y[0] = YC[i]                                                                # Panel starting Y point
    Y[1] = YC[i] + S[i]*np.sin(delta[i])                                        # Panel ending Y point
    if (i == 0):                                                                # For first panel
        plt.plot(X,Y,'b-',label='First Panel')                                  # Plot the first panel normal vector
    elif (i == 1):                                                              # For second panel
        plt.plot(X,Y,'g-',label='Second Panel')                                 # Plot the second panel normal vector
    else:                                                                       # For every other panel
        plt.plot(X,Y,'r-')                                                      # Plot the panel normal vector
plt.xlabel('X-Axis')                                                            # Set X-label
plt.ylabel('Y-Axis')                                                            # Set Y-label
plt.title('Panel Geometry')                                                     # Set title
plt.axis('equal')                                                               # Set axes equal
plt.legend()                                                                    # Plot legend
plt.show()                                                                      # Display plot
