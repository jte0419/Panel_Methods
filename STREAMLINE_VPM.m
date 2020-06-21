function [Nx,Ny] = STREAMLINE_VPM(XP,YP,XB,YB,phi,S)

% FUNCTION - COMPUTE Nx AND Ny GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 01/23/19
% Updated   : 01/23/19 - Started code
%                      - Works as expected
%           : 04/28/20 - Updated E value error handling to match Python
% 
% PURPOSE
% - Compute the integral expression for constant strength vortex panels
% - Vortex panel strengths are constant, but can change from panel to panel
% - Geometric integral for X-velocity: Nx(pj)
% - Geometric integral for Y-velocity: Ny(pj)
% 
% REFERENCES
% - [1]: Streamline Geometric Integral VPM, Nx(pj) and Ny(pj)
%           Link: https://www.youtube.com/watch?v=TBwBnW87hso
% INPUTS
% - XP     : X-coordinate of computation point, P
% - YP     : Y-coordinate of computation point, P
% - XB     : X-coordinate of boundary points
% - YB     : Y-coordinate of boundary points
% - phi    : Angle between positive X-axis and interior of panel
% - S      : Length of panel
% 
% OUTPUTS
% - Nx     : Value of X-direction geometric integral
% - Ny     : Value of Y-direction geometric integral

% Number of panels
numPan = length(XB)-1;                                                      % Number of panels (control points)

% Initialize arrays
Nx = zeros(numPan,1);                                                       % Initialize Nx integral array
Ny = zeros(numPan,1);                                                       % Initialize Ny integral array

% Compute Nx and Ny
for j = 1:1:numPan                                                          % Loop over all panels
    % Compute intermediate values
    A  = -(XP-XB(j))*cos(phi(j))-(YP-YB(j))*sin(phi(j));                    % A term
    B  = (XP-XB(j))^2+(YP-YB(j))^2;                                         % B term
    Cx = sin(phi(j));                                                       % Cx term (X-direction)
    Dx = -(YP-YB(j));                                                       % Dx term (X-direction)
    Cy = -cos(phi(j));                                                      % Cy term (Y-direction)
    Dy = XP-XB(j);                                                          % Dy term (Y-direction)
    E  = sqrt(B-A^2);                                                       % E term
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Nx
    term1 = 0.5*Cx*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Nx equation
    term2 = ((Dx-A*Cx)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Nx equation
    Nx(j) = term1 + term2;                                                  % Compute Nx integral
    
    % Compute Ny
    term1 = 0.5*Cy*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Ny equation
    term2 = ((Dy-A*Cy)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Ny equation
    Ny(j) = term1 + term2;                                                  % Compute Ny integral
    
	% Zero out any NANs, INFs, or imaginary numbers
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
end
