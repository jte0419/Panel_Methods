function [K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S)

% FUNCTION - COMPUTE K AND L GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 01/23/19
% Updated   : 01/23/19 - Started code
%                      - Works as expected
%           : 04/28/20 - Updated E value error handling to match Python
% 
% PUROSE
% - Compute the integral expression for constant strength vortex panels
% - Vortex panel strengths are constant, but can change from panel to panel
% - Geometric integral for panel-normal    : K(ij)
% - Geometric integral for panel-tangential: L(ij)
%
% REFERENCES
% - [1]: Normal Geometric Integral VPM, K(ij)
%           Link: https://www.youtube.com/watch?v=5lmIv2CUpoc
% - [2]: Tangential Geometric Integral VPM, L(ij)
%           Link: https://www.youtube.com/watch?v=IxWJzwIG_gY
%
% INPUTS
% - XC  : X-coordinate of control points
% - YC  : Y-coordinate of control points
% - XB  : X-coordinate of boundary points
% - YB  : Y-coordinate of boundary points
% - phi : Angle between positive X-axis and interior of panel
% - S   : Length of panel
% 
% OUTPUTS
% - K   : Value of panel-normal integral (Ref [1])
% - L   : Value of panel-tangential integral (Ref [2])

% Number of panels
numPan = length(XC);                                                        % Number of panels

% Initialize arrays
K = zeros(numPan,numPan);                                                   % Initialize K integral matrix
L = zeros(numPan,numPan);                                                   % Initialize L integral matrix

% Compute integral
for i = 1:1:numPan                                                          % Loop over i panels
    for j = 1:1:numPan                                                      % Loop over j panels
        if (j ~= i)                                                         % If panel j is not the same as panel i
            A  = -(XC(i)-XB(j))*cos(phi(j))-(YC(i)-YB(j))*sin(phi(j));      % A term
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;                           % B term
            Cn = -cos(phi(i)-phi(j));                                       % C term (normal)
            Dn = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j))*sin(phi(i));       % D term (normal)
            Ct = sin(phi(j)-phi(i));                                        % C term (tangential)
            Dt = (XC(i)-XB(j))*sin(phi(i))-(YC(i)-YB(j))*cos(phi(i));       % D term (tangential)
            E  = sqrt(B-A^2);                                               % E term
            if (~isreal(E))
                E = 0;
            end
            
            % Compute K
            term1  = 0.5*Cn*log((S(j)^2+2*A*S(j)+B)/B);                     % First term in K equation
            term2  = ((Dn-A*Cn)/E)*(atan2((S(j)+A),E)-atan2(A,E));          % Second term in K equation
            K(i,j) = term1 + term2;                                         % Compute K integral
            
            % Compute L
            term1  = 0.5*Ct*log((S(j)^2+2*A*S(j)+B)/B);                     % First term in L equation
            term2  = ((Dt-A*Ct)/E)*(atan2((S(j)+A),E)-atan2(A,E));          % Second term in L equation
            L(i,j) = term1 + term2;                                         % Compute L integral
        end
        
        % Zero out any NANs, INFs, or imaginary numbers
        if (isnan(K(i,j)) || isinf(K(i,j)) || ~isreal(K(i,j)))
            K(i,j) = 0;
        end
        if (isnan(L(i,j)) || isinf(L(i,j)) || ~isreal(L(i,j)))
            L(i,j) = 0;
        end
    end
end
