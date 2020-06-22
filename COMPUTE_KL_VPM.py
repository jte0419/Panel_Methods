# FUNCTION - COMPUTE K AND L GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD
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
# - Geometric integral for panel-normal    : K(ij)
# - Geometric integral for panel-tangential: L(ij)
#
# REFERENCES
# - [1]: Normal Geometric Integral VPM, K(ij)
#           Link: https://www.youtube.com/watch?v=5lmIv2CUpoc
# - [2]: Tangential Geometric Integral VPM, L(ij)
#           Link: https://www.youtube.com/watch?v=IxWJzwIG_gY
#
# INPUTS
# - XC  : X-coordinate of control points
# - YC  : Y-coordinate of control points
# - XB  : X-coordinate of boundary points
# - YB  : Y-coordinate of boundary points
# - phi : Angle between positive X-axis and interior of panel
# - S   : Length of panel
# 
# OUTPUTS
# - K   : Value of panel-normal integral (Ref [1])
# - L   : Value of panel-tangential integral (Ref [2])

import numpy as np
import math as math

np.seterr('raise')

def COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S):
    
    # Number of panels
    numPan = len(XC)                                                                # Number of panels
    
    # Initialize arrays
    K = np.zeros([numPan,numPan])                                                   # Initialize K integral matrix
    L = np.zeros([numPan,numPan])                                                   # Initialize L integral matrix
    
    # Compute integral
    for i in range(numPan):                                                         # Loop over i panels
        for j in range(numPan):                                                     # Loop over j panels
            if (j != i):                                                            # If panel j is not the same as panel i
                # Compute intermediate values
                A  = -(XC[i]-XB[j])*np.cos(phi[j])-(YC[i]-YB[j])*np.sin(phi[j])     # A term
                B  = (XC[i]-XB[j])**2 + (YC[i]-YB[j])**2                            # B term
                Cn = -np.cos(phi[i]-phi[j])                                         # C term (normal)
                Dn = (XC[i]-XB[j])*np.cos(phi[i])+(YC[i]-YB[j])*np.sin(phi[i])      # D term (normal)
                Ct = np.sin(phi[j]-phi[i])                                          # C term (tangential)
                Dt = (XC[i]-XB[j])*np.sin(phi[i])-(YC[i]-YB[j])*np.cos(phi[i])      # D term (tangential)
                E  = np.sqrt(B-A**2)                                                # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):       # If E term is 0 or complex or a NAN or an INF
                    K[i,j] = 0                                                      # Set K value equal to zero
                    L[i,j] = 0                                                      # Set L value equal to zero
                else:
                    # Compute K
                    term1  = 0.5*Cn*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in K equation
                    term2  = ((Dn-A*Cn)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in K equation
                    K[i,j] = term1 + term2                                          # Compute K integral
                    
                    # Compute L
                    term1  = 0.5*Ct*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in L equation
                    term2  = ((Dt-A*Ct)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in L equation
                    L[i,j] = term1 + term2                                          # Compute L integral
            
            # Zero out any problem values
            if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):      # If K term is complex or a NAN or an INF
                K[i,j] = 0                                                          # Set K value equal to zero
            if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):      # If L term is complex or a NAN or an INF
                L[i,j] = 0                                                          # Set L value equal to zero
    
    return K, L                                                                     # Return both K and L matrices
