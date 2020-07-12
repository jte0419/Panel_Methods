function [Mx,My] = STREAMLINE_SPM_N(XP,YP,XB,YB,phi,S,numPan,jInd)

% FUNCTION - COMPUTE MX AND MY GEOMETRIC INTEGRALS FOR SOURCE PANEL METHOD // N AIRFOILS
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% 
% PURPOSE
% - Compute the geometric integral at point P due to source panels
% - Source panel strengths are constant, but can change from panel to panel
% - Geometric integral for X-direction: Mx(pj)
% - Geometric integral for Y-direction: My(pj)
% 
% REFERENCE
% - [1]: Streamline Geometric Integral SPM, Mx(pj) and My(pj)
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
% - Mx     : Value of X-direction geometric integral (Ref [1])
% - My     : Value of Y-direction geometric integral (Ref [1])

% Initialize arrays
Mx = zeros(numPan,1);                                                       % Initialize Mx integral array
My = zeros(numPan,1);                                                       % Initialize My integral array

% Compute Mx and My
for j = 1:1:numPan                                                          % Loop over the j panels
    % Compute intermediate values
    A  = -(XP-XB(jInd(j)))*cos(phi(jInd(j))) - ...                          % A term
          (YP-YB(jInd(j)))*sin(phi(jInd(j)));
    B  = (XP-XB(jInd(j)))^2+(YP-YB(jInd(j)))^2;                             % B term
    Cx = -cos(phi(jInd(j)));                                                % C term (X-direction)
    Dx = XP - XB(jInd(j));                                                  % D term (X-direction)
    Cy = -sin(phi(jInd(j)));                                                % C term (Y-direction)
    Dy = YP - YB(jInd(j));                                                  % D term (Y-direction)
    E  = sqrt(B-A^2);                                                       % E term
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Mx, Ref [1]
    term1 = 0.5*Cx*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);                  % First term in Mx equation
    term2 = ((Dx-A*Cx)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));           % Second term in Mx equation
    Mx(j) = term1 + term2;                                                  % X-direction geometric integral
    
    % Compute My, Ref [1]
    term1 = 0.5*Cy*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);                  % First term in My equation
    term2 = ((Dy-A*Cy)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));           % Second term in My equation
    My(j) = term1 + term2;                                                  % Y-direction geometric integral
    
    % Zero out any NANs, INFs, or imaginary numbers
    if (isnan(Mx(j)) || isinf(Mx(j)) || ~isreal(Mx(j)))
        Mx(j) = 0;
    end
    if (isnan(My(j)) || isinf(My(j)) || ~isreal(My(j)))
        My(j) = 0;
    end
end
