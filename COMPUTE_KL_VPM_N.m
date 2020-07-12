function [K,L] = COMPUTE_KL_VPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd)

% FUNCTION - COMPUTE K AND L GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD // N AIRFOILS
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 01/23/19
% Updated   : 01/23/19 - Started code
%                      - Works as expected
%             04/16/19 - Adding N-airfoil capability
% 
% PUROSE
% - Compute the integral expression for vortex panel method
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
% - XC     : X-coordinate of control points
% - YC     : Y-coordinate of control points
% - XB     : X-coordinate of boundary points
% - YB     : Y-coordinate of boundary points
% - phi    : Angle between positive X-axis and interior of panel
% - S      : Length of panel
% - numPan : Number of panels
% - iInd   : Actual i panel indices (not inter-airfoil panels)
% - jInd   : Actual j panel indices (not inter-airfoil panels)
% 
% OUTPUTS
% - K   : Value of panel-normal integral (Ref [1])
% - L   : Value of panel-tangential integral (Ref [2])

% Initialize arrays
K = zeros(numPan,numPan);                                                   % Initialize K integral matrix
L = zeros(numPan,numPan);                                                   % Initialize L integral matrix

% Compute integral
for i = 1:1:numPan                                                          % Loop over i panels
    for j = 1:1:numPan                                                      % Loop over j panels
        if (jInd(j) ~= iInd(i))                                             % If panel j is not the same as panel i
            A  = -(XC(iInd(i))-XB(jInd(j)))*cos(phi(j)) - ...               % A term
                  (YC(iInd(i))-YB(jInd(j)))*sin(phi(j));
            B  = (XC(iInd(i))-XB(jInd(j)))^2+(YC(iInd(i))-YB(jInd(j)))^2;   % B term
            Cn = -cos(phi(iInd(i))-phi(jInd(j)));                           % C term (normal)
            Dn = (XC(iInd(i))-XB(jInd(j)))*cos(phi(iInd(i))) + ...          % D term (normal)
                 (YC(iInd(i))-YB(jInd(j)))*sin(phi(iInd(i)));
            Ct = sin(phi(jInd(j))-phi(iInd(i)));                            % C term (tangential)
            Dt = (XC(iInd(i))-XB(jInd(j)))*sin(phi(iInd(i))) - ...          % D term (tangential)
                 (YC(iInd(i))-YB(jInd(j)))*cos(phi(iInd(i)));
            E  = sqrt(B-A^2);                                               % E term
            if (isnan(E) || ~isreal(E))                                     % If E is a NaN or not real
                E = 0;                                                      % Set E equal to zero
            end
            
            % Compute K
            term1 = 0.5*Cn*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);          % First term in K equation
            term2 = ((Dn-A*Cn)/E)*(atan2((S(jInd(j))+A),E)-atan2(A,E));     % Second term in K equation
            K(i,j) = term1 + term2;                                         % Compute K integral
            
            % Compute L
            term1 = 0.5*Ct*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);          % First term in L equation
            term2 = ((Dt-A*Ct)/E)*(atan2((S(jInd(j))+A),E)-atan2(A,E));     % Second term in L equation
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
