# FUNCTION - COMPUTE CIRCULATION
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started: 02/19/19
# Updated: 02/19/19 - Transferred from MATLAB to Python
#                   - Works as expected
#
# PURPOSE
# - Compute the circulation around the defined ellipse
# 
# INPUTS
# - a    : Horizontal axis half-length
# - b    : Vertical axis half-length
# - x0   : Ellipse center X coordinate
# - y0   : Ellipse center Y coordinate
# - numT : Number of points for integral
# - XX   : Meshgrid X values
# - YY   : Meshgrid Y values
#
# OUTPUTS
# - Gamma : Circulation [length^2/time]
# - xC    : X-values of integral curve [numT x 1]
# - yC    : Y-values of integral curve [numT x 1]
# - VxC   : Velocity X-component on integral curve [numT x 1]
# - VyC   : Velocity Y-component on integral curve [numT x 1]

import numpy as np
from scipy import interpolate

def COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,X,Y):
    
    t     = np.linspace(0,2*np.pi,numT)                                         # Discretized ellipse into angles [rad]
    xC    = a*np.cos(t) + x0                                                    # X coordinates of ellipse
    yC    = b*np.sin(t) + y0                                                    # Y coordinates of ellipse
    fx    = interpolate.RectBivariateSpline(Y,X,Vx)                             # Interpolate X velocities from grid to ellipse points
    fy    = interpolate.RectBivariateSpline(Y,X,Vy)                             # Interpolate Y velocities from grid to ellipse points
    VxC   = fx.ev(yC,xC)                                                        # X velocity component on ellipse
    VyC   = fy.ev(yC,xC)                                                        # Y velocity component on ellipse
    Gamma = -(np.trapz(VxC,xC) + np.trapz(VyC,yC))                              # Compute integral using trapezoid rule
    
    return Gamma, xC, yC, VxC, VyC
