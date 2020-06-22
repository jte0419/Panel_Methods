# FUNCTION - COMPUTE Nx AND Ny GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 01/23/19
# Updated   : 01/23/19 - Started code in MATLAB
#                      - Works as expected
#           : 02/03/19 - Transferred to Python
#                      - Works as expected
#           : 04/28/20 - Fixed E value error handling
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

np.seterr('raise')

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
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):           # If E term is 0 or complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
            Ny[j] = 0                                                           # Set Ny value equal to zero
        else:
            # Compute Nx, Ref [1]
            term1 = 0.5*Cx*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Nx equation
            term2 = ((Dx-A*Cx)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Nx equation
            Nx[j] = term1 + term2;                                              # Compute Nx integral
            
            # Compute Ny, Ref [1]
            term1 = 0.5*Cy*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Ny equation
            term2 = ((Dy-A*Cy)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Ny equation
            Ny[j] = term1 + term2;                                              # Compute Ny integral
            
        # Zero out any problem values
        if (np.iscomplex(Nx[j]) or np.isnan(Nx[j]) or np.isinf(Nx[j])):         # If Nx term is complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
        if (np.iscomplex(Ny[j]) or np.isnan(Ny[j]) or np.isinf(Ny[j])):         # If Ny term is complex or a NAN or an INF
            Ny[j] = 0                                                           # Set Ny value equal to zero
    
    return Nx, Ny                                                               # Return both Nx and Ny matrices
