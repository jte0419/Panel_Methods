# CALLING XFOIL FROM PYTHON
# Written by: JoshTheEngineer
# YouTube: www.youtube.com/joshtheengineer
# Website: www.joshtheengineer.com
# Started: 01/01/19
# Updated: 01/01/19 - Started code in MATLAB
#                   - Works as expected
#          02/03/19 - Transferred code from MATLAB to Python
#                   - Works as expected

import os
import numpy as np
import matplotlib.pyplot as plt

# %% CREATE LOADING FILE

# Knowns
NACA       = '0012'
AoA        = '0'
numNodes   = '170'
saveFlnmAF = 'Save_Airfoil.txt'
saveFlnmCp = 'Save_Cp.txt'
xfoilFlnm  = 'xfoil_input.txt'

# Delete files if they exist
if os.path.exists(saveFlnmAF):
    os.remove(saveFlnmAF)

if os.path.exists(saveFlnmCp):
    os.remove(saveFlnmCp)
    
# Create the airfoil
fid = open(xfoilFlnm,"w")
fid.write("NACA " + NACA + "\n")
fid.write("PPAR\n")
fid.write("N " + numNodes + "\n")
fid.write("\n\n")
fid.write("PSAV " + saveFlnmAF + "\n")
fid.write("OPER\n")
fid.write("ALFA " + AoA + "\n")
fid.write("CPWR " + saveFlnmCp + "\n")
fid.close()

# Run the XFoil calling command
os.system("xfoil.exe < xfoil_input.txt")

# Delete file after running
if os.path.exists(xfoilFlnm):
    os.remove(xfoilFlnm)

# %% READ DATA FILE: AIRFOIL

flpth = "C:/Users/Josh/Documents/Python/Panel_Methods/"
flnm  = flpth + saveFlnmAF
    
# Load the data from the text file
dataBuffer = np.loadtxt(flnm, skiprows=0)

# Extract data from the loaded dataBuffer array
XB = dataBuffer[:,0]
YB = dataBuffer[:,1]

# Delete file after loading
if os.path.exists(saveFlnmAF):
    os.remove(saveFlnmAF)

# %% READ DATA FILE: PRESSURE COEFFICIENT

# Load the data from the text file
dataBuffer = np.loadtxt(saveFlnmCp, skiprows=3)

# Extract data from the loaded dataBuffer array
X_0  = dataBuffer[:,0]
Y_0  = dataBuffer[:,1]
Cp_0 = dataBuffer[:,2]

# Delete file after loading
if os.path.exists(saveFlnmCp):
    os.remove(saveFlnmCp)

# %% EXTRACT UPPER AND LOWER AIRFOIL DATA

# Split airfoil into (U)pper and (L)ower
XB_U = XB[YB >= 0]
XB_L = XB[YB < 0]
YB_U = YB[YB >= 0]
YB_L = YB[YB < 0]

# Split XFoil results into (U)pper and (L)ower
Cp_U = Cp_0[YB >= 0]
Cp_L = Cp_0[YB < 0]
X_U  = X_0[YB >= 0]
X_L  = X_0[YB < 0]

# %% PLOT DATA

# Plot airfoil
fig = plt.figure(1)
plt.cla()
plt.plot(XB_U,YB_U,'b.-',label='Upper')
plt.plot(XB_L,YB_L,'r.-',label='Lower')
plt.xlabel('X-Coordinate')
plt.ylabel('Y-Coordinate')
plt.title('Airfoil')
plt.axis('equal')
plt.legend()
plt.show()

# Plot pressure coefficient
fig = plt.figure(2)
plt.cla()
plt.plot(X_U,Cp_U,'b.-',label='Upper')
plt.plot(X_L,Cp_L,'r.-',label='Lower')
plt.xlim(0,1)
plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.title('Pressure Coefficient')
plt.show()
plt.legend()
plt.gca().invert_yaxis()
