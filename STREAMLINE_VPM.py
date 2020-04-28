# FUNCTION - COMPUTE Nx AND Ny GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD
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
# - Compute the integral expression for constant strength vortex panels
# - Vortex panel strengths are constant, but can change from panel to panel
# - Geometric integral for X-direction: Nx(pj)
# - Geometric integral for Y-direction: Ny(pj)
#
# REFERENCES
# - [1]: Streamline Geometric Integral VPM, Nx(pj) and Ny(pj)
#           Link: https://www.youtube.com/watch?v=TBwBnW87hso
#
# INPUTS
# - XP  : X-coordinate of computation point, P
# - YP  : Y-coordinate of computation point, P
# - XB  : X-coordinate of boundary points
# - YB  : Y-coordinate of boundary points
# - phi : Angle between positive X-axis and interior of panel
# - S   : Length of panel
# 
# OUTPUTS
# - Nx  : Value of X-direction geometric integral
# - Ny  : Value of Y-direction geometric integral

import numpy as np
import math as math

def STREAMLINE_VPM(XP,YP,XB,YB,phi,S):
    
    # Number of panels
    numPan = len(XB)-1                                                          # Number of panels (control points)
    
    # Initialize arrays
    Nx = np.zeros(numPan)                                                       # Initialize Nx integral array
    Ny = np.zeros(numPan)                                                       # Initialize Ny integral array
    
    # Compute Nx and Ny
    for j in range(numPan):                                                     # Loop over all panels
        # Compute intermediate values
        A = -(XP-XB[j])*np.cos(phi[j]) - (YP-YB[j])*np.sin(phi[j])              # A term
        B  = (XP-XB[j])**2 + (YP-YB[j])**2                                      # B term
        Cx = np.sin(phi[j])                                                     # Cx term (X-direction)
        Dx = -(YP-YB[j])                                                        # Dx term (X-direction)
        Cy = -np.cos(phi[j])                                                    # Cy term (Y-direction)
        Dy = XP-XB[j]                                                           # Dy term (Y-direction)
        E  = math.sqrt(B-A**2)                                                  # E term
        if (np.isnan(E) or ~np.isreal(E)):                                      # If E is a NaN or imaginary
            E = 0                                                               # Set E equal to zero
                    
        # Compute Nx
        term1 = 0.5*Cx*np.log((S[j]**2 + 2*A*S[j] + B)/B);                      # First term in Nx equation
        term2 = ((Dx-A*Cx)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));       # Second term in Nx equation
        Nx[j] = term1 + term2;                                                  # Compute Nx integral
        
        # Compute Ny
        term1 = 0.5*Cy*np.log((S[j]**2 + 2*A*S[j]+B)/B);                        # First term in Ny equation
        term2 = ((Dy-A*Cy)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));       # Second term in Ny equation
        Ny[j] = term1 + term2;                                                  # Compute Ny integral
    
    return Nx, Ny                                                               # Return both Nx and Ny matrices