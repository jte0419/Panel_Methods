# FUNCTION - LOAD AIRFOIL (SELIG FORMAT)
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 01/23/19
# Updated   : 01/23/19 - Started code in MATLAB
#                      - Works as expected
#           : 02/03/19 - Transferred to Python
#                      - Works as expected
#
# PURPOSE
# - Load the X and Y coordinates for airfoil from fileName
# - Data comes from files on UIUC website
# - Loads every airfoil except naca23015.dat (has weird data)
# - Checked to make sure the loaded data works with panel method codes
#
# INPUTS
# - fileName : File name of the airfoil without extension (must be .dat)
# 
# OUTPUTS
# - dataX : Airfoil X-coordinate data array [Nx1]
# - dataY : Airfoil Y-coordinate data array [Nx2]

import numpy as np

def LOAD_AIRFOIL_SELIG(fileName):
    
    # Define the number of header lines in the file
    hdrlns = 1
    if (fileName == 'nasasc2-0714'):
        hdrlns = 3
    elif (fileName == 's1020'):
        hdrlns = 2
    
    # Load the data from the text file
    flpth = "C:/Users/Josh/Documents/MATLAB/Panel_Methods/Panel_Method_Complete/Airfoil_DAT_Selig/"
    flnm = flpth + fileName + ".dat"
    dataBuffer = np.loadtxt(flnm,delimiter=' ', skiprows=hdrlns)
    
    # Extract data from the loaded dataBuffer array
    dataX = dataBuffer[:,0]
    dataY = dataBuffer[:,1]
    
    return dataX, dataY