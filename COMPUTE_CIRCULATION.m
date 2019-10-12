% FUNCTION: COMPUTE CIRCULATION
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 02/11/19
% Updated: 02/11/19 - Started code
%                   - Works as expected
%
% PURPOSE
% - Compute the circulation around the defined ellipse
% 
% INPUTS
% - a    : Horizontal axis half-length
% - b    : Vertical axis half-length
% - x0   : Ellipse center X coordinate
% - y0   : Ellipse center Y coordinate
% - numT : Number of points for integral
% - XX   : Meshgrid X values
% - YY   : Meshgrid Y values
%
% OUTPUTS
% - Lambda : Circulation [length^2/time]
% - xC     : X-values of integral curve [numT x 1]
% - yC     : Y-values of integral curve [numT x 1]
% - VxC    : Velocity X-component on integral curve [numT x 1]
% - VyC    : Velocity Y-component on integral curve [numT x 1]

function [Gamma,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY)

tEnd  = (2*pi) - ((2*pi)/numT);                                             % Ending angle [rad]
t     = linspace(0,tEnd,numT)';                                             % Discretize ellipse into angles [rad]
xC    = a*cos(t) + x0;                                                      % X coordinates of ellipse
yC    = b*sin(t) + y0;                                                      % Y coordinates of ellipse
VxC   = interp2(XX,YY,Vx,xC,yC);                                            % X velocity component on ellipse
VyC   = interp2(XX,YY,Vy,xC,yC);                                            % Y velocity component on ellipse
Gamma = -(trapz(xC,VxC) + trapz(yC,VyC));                                   % Compute integral using trapezoid rule

