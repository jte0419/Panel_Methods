function [I,J] = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S)

% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Updated   : 04/28/20 - Updated E value error handling to match Python
%
% PURPOSE
% - Compute the integral expression for constant strength source panels
% - Source panel strengths are constant, but can change from panel to panel
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
% - XC  : X-coordinate of control points
% - YC  : Y-coordinate of control points
% - XB  : X-coordinate of boundary points
% - YB  : Y-coordinate of boundary points
% - phi : Angle between positive X-axis and interior of panel
% - S   : Length of panel
% 
% OUTPUTS
% - I   : Value of panel-normal integral (Eq. 3.163 in Anderson or Ref [1])
% - J   : Value of panel-tangential integral (Eq. 3.165 in Anderson or Ref [2])

% Number of panels
numPan = length(XC);                                                        % Number of panels/control points

% Initialize arrays
I = zeros(numPan,numPan);                                                   % Initialize I integral matrix
J = zeros(numPan,numPan);                                                   % Initialize J integral matrix

% Compute integral
for i = 1:1:numPan                                                          % Loop over i panels
    for j = 1:1:numPan                                                      % Loop over j panels
        if (j ~= i)                                                         % If the i and j panels are not the same
            % Compute intermediate values
            A  = -(XC(i)-XB(j))*cos(phi(j))-(YC(i)-YB(j))*sin(phi(j));      % A term
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;                           % B term
            Cn = sin(phi(i)-phi(j));                                        % C term (normal)
            Dn = -(XC(i)-XB(j))*sin(phi(i))+(YC(i)-YB(j))*cos(phi(i));      % D term (normal)
            Ct = -cos(phi(i)-phi(j));                                       % C term (tangential)
            Dt = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j))*sin(phi(i));       % D term (tangential)
            E  = sqrt(B-A^2);                                               % E term
            if (~isreal(E))
                E = 0;
            end
            
            % Compute I (needed for normal velocity), Ref [1]
            term1  = 0.5*Cn*log((S(j)^2+2*A*S(j)+B)/B);                     % First term in I equation
            term2  = ((Dn-A*Cn)/E)*(atan2((S(j)+A),E) - atan2(A,E));        % Second term in I equation
            I(i,j) = term1 + term2;                                         % Compute I integral
            
            % Compute J (needed for tangential velocity), Ref [2]
            term1  = 0.5*Ct*log((S(j)^2+2*A*S(j)+B)/B);                     % First term in J equation
            term2  = ((Dt-A*Ct)/E)*(atan2((S(j)+A),E) - atan2(A,E));        % Second term in J equation
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
