function [Nx,Ny] = STREAMLINE_VPM_N(XP,YP,XB,YB,phi,S,numPan,jInd)

% FUNCTION - COMPUTE NX AND NY GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD // N AIRFOILS
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% 
% PURPOSE
% - Compute the geometric integral at point P due to source panels
% - Source panel strengths are constant, but can change from panel to panel
% - Geometric integral for X-direction: Nx(pj)
% - Geometric integral for Y-direction: Ny(pj)
% 
% REFERENCE
% - [1]: Streamline Geometric Integral SPM, Nx(pj) and Ny(pj)
%           Link: https://www.youtube.com/watch?v=BnPZjGCatcg
% 
% INPUTS
% - XP     : X-coordinate of computation point, P
% - YP     : Y-coordinate of computation point, P
% - XB     : X-coordinate of boundary points
% - YB     : Y-coordinate of boundary points
% - phi    : Angle between positive X-axis and interior of panel
% - S      : Length of panel
% - numPan : Number of panels
% - jInd   : Actual panel indices (not inter-airfoil panels)
% 
% OUTPUTS
% - Nx     : Value of X-direction geometric integral (Ref [1])
% - Ny     : Value of Y-direction geometric integral (Ref [1])

% Initialize arrays
Nx = zeros(numPan,1);                                                       % Initialize Nx integral array
Ny = zeros(numPan,1);                                                       % Initialize Ny integral array

% Compute Nx and Ny
for j = 1:1:numPan                                                          % Loop over all panels
    % Compute intermediate values
    A  = -(XP-XB(jInd(j)))*cos(phi(jInd(j))) - ...                          % A term
          (YP-YB(jInd(j)))*sin(phi(jInd(j)));  
    B  = (XP-XB(jInd(j)))^2+(YP-YB(jInd(j)))^2;                             % B term
    Cx = sin(phi(jInd(j)));                                                 % Cx term (X-direction)
    Dx = -(YP-YB(jInd(j)));                                                 % Dx term (X-direction)
    Cy = -cos(phi(jInd(j)));                                                % Cy term (Y-direction)
    Dy = XP-XB(jInd(j));                                                    % Dy term (Y-direction)
    E  = sqrt(B-A^2);                                                       % E term
    if (~isreal(E))
        E = 0;
    end
	
    % Compute Nx
    term1 = 0.5*Cx*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);                  % First term in Nx equation
    term2 = ((Dx-A*Cx)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));           % Second term in Nx equation
    Nx(j) = term1 + term2;                                                  % Compute Nx integral
    
    % Compute Ny
    term1 = 0.5*Cy*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);                  % First term in Ny equation
    term2 = ((Dy-A*Cy)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));           % Second term in Ny equation
    Ny(j) = term1 + term2;                                                  % Compute Ny integral
    
	% Zero out any NANs, INFs, or imaginary numbers
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
end

