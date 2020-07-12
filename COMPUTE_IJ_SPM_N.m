function [I,J] = COMPUTE_IJ_SPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd)

% FUNCTION - COMPUTE I AND J GEOMETRIC INTEGRALS FOR SOURCE PANEL METHOD // N AIRFOILS
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Updated   : 04/16/19 - Adding N-airfoil capability
% 
% PURPOSE
% - Compute the integral expression for source panel method
% - Geometric integral for panel-normal    : I(ij)
% - Geometric integral for panel-tangential: J(ij)
% 
% REFERENCES
% - [1]: Normal Geometric Integral SPM, I(ij)
%           Link: https://www.youtube.com/watch?v=76vPudNET6U
% - [2]: Tangential Geometric Integral SPM, J(ij)
%           Link: https://www.youtube.com/watch?v=JRHnOsueic8
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
% - I   : Value of panel-normal integral (Eq. 3.163 in Anderson or Ref [1])
% - J   : Value of panel-tangential integral (Eq. 3.165 in Anderson or Ref [2])

% Initialize arrays
I = zeros(numPan,numPan);                                                   % Initialize I integral matrix
J = zeros(numPan,numPan);                                                   % Initialize J integral matrix

% Compute integral
for i = 1:1:numPan                                                          % Loop over i panels
    for j = 1:1:numPan                                                      % Loop over j panels
        if (jInd(j) ~= iInd(i))                                             % If the i and j panels are not the same
            % Compute intermediate values
            A  = -(XC(iInd(i))-XB(jInd(j)))*cos(phi(jInd(j))) - ...         % A term
                  (YC(iInd(i))-YB(jInd(j)))*sin(phi(jInd(j)));      
            B  = (XC(iInd(i))-XB(jInd(j)))^2+(YC(iInd(i))-YB(jInd(j)))^2;   % B term
            Cn = sin(phi(iInd(i))-phi(jInd(j)));                            % C term (normal)
            Dn = -(XC(iInd(i))-XB(jInd(j)))*sin(phi(iInd(i))) + ...         % D term (normal)
                  (YC(iInd(i))-YB(jInd(j)))*cos(phi(iInd(i)));
            Ct = -cos(phi(iInd(i))-phi(jInd(j)));                           % C term (tangential)
            Dt = (XC(iInd(i))-XB(jInd(j)))*cos(phi(iInd(i))) + ...          % D term (tangential)
                 (YC(iInd(i))-YB(jInd(j)))*sin(phi(iInd(i)));
            E  = sqrt(B-A^2);                                               % E term
            if (~isreal(E))                                                 % If E is imaginary
                E = 0;                                                      % Set E equal to zero
            end
            
            % Compute I (needed for normal velocity), Ref [1]
            term1  = 0.5*Cn*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);         % First term in I equation
            term2  = ((Dn-A*Cn)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));  % Second term in I equation
            I(i,j) = term1 + term2;                                         % Compute I integral
            
            % Compute J (needed for tangential velocity), Ref [2]
            term1  = 0.5*Ct*log((S(jInd(j))^2+2*A*S(jInd(j))+B)/B);         % First term in J equation
            term2  = ((Dt-A*Ct)/E)*(atan2((S(jInd(j))+A),E) - atan2(A,E));  % Second term in J equation
            J(i,j) = term1 + term2;                                         % Compute J integral
        end
        
        % Zero out any NANs, INFs, or imaginary numbers
        if (isnan(I(i,j)) || isinf(I(i,j)) || ~isreal(I(i,j)))
            I(i,j) = 0;
        end
        if (isnan(J(i,j)) || isinf(J(i,j)) || ~isreal(J(i,j)))
            J(i,j) = 0;
        end
    end
end
