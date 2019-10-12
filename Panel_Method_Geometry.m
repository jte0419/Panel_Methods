% PANEL METHOD GEOMETRY
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 01/24/19
% Updated: 01/24/19 - Started code
%                   - Works as expected

clear;
clc;

%% CREATE/LOAD GEOMETRY

% Load circle or airfoil
% load = 'Airfoil';                                                           % Use the airfoil geometry
load = 'Circle';                                                            % Use the circle geometry
AoA  = 0;                                                                   % Angle of attack [deg]

if (strcmpi(load,'Airfoil'))                                                % If airfoil is to be loaded
    fileName = 'goe623';                                                    % Specify the airfoil filename
    [XB,YB]  = LOAD_AIRFOIL_SELIG(fileName);                                % Load X and Y data points of airfoil
elseif (strcmpi(load,'Circle'))                                             % If circle is to be used
    numB = 9;                                                               % Number of boundary points (including endpt)
    tO   = 22.5;                                                            % Boundary point angle offset [deg]
    
    % Angles used to compute boundary points
    theta = linspace(0,360,numB)';                                          % Create angles for computing boundary point locations [deg]
    theta = theta + tO;                                                     % Add angle offset [deg]
    theta = theta*(pi/180);                                                 % Convert from degrees to radians [rad]
    
    % Boundary points
    XB = cos(theta);                                                        % Compute boundary point X-value (radius of 1)
    YB = sin(theta);                                                        % Compute boundary point Y-value (radius of 1)
end

% Number of panels
numPan = length(XB)-1;                                                      % Number of panels

% Check for direction of points
edge = zeros(numPan,1);                                                     % Initialize edge value array
for i = 1:1:numPan                                                          % Loop over all panels
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              % Compute edge value
end
sumEdge = sum(edge);                                                        % Sum of all edge values

% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            % If panels are CCW
    fprintf('Points are counter-clockwise.  Flipping.\n');                  % Display message
    XB = flipud(XB);                                                        % Flip the X-data array
    YB = flipud(YB);                                                        % Flip the Y-data array
elseif (sumEdge > 0)                                                        % If the panels are CW already
    fprintf('Points are clockwise.  Not flipping.\n');                      % Display message
end

%% COMPUTE GEOMETRIC VARIABLES

% Find geometric quantities of airfoil
XC   = zeros(numPan,1);                                                     % Initialize control point X-coordinate array
YC   = zeros(numPan,1);                                                     % Initialize control point Y-coordinate array
S    = zeros(numPan,1);                                                     % Initialize panel length array
phiD = zeros(numPan,1);                                                     % Initialize panel orientation angle array [deg]
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (right horizontal to inside face)
    if (phiD(i) < 0)                                                        % Make all angles positive
        phiD(i) = phiD(i) + 360;                                            % Add 360 degrees to angle
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle of panel normal
betaD              = deltaD - AoA;                                          % Angle of panel normal and AoA [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make sure angles aren't greater than 360 [deg]

%% PLOTTING

T = linspace(0,360,1000)';                                                  % Angle array to compute dashed circle
x = cosd(T);                                                                % Circle X points
y = sind(T);                                                                % Circle Y points

% Plot shape polygon with panel normal vectors
figure(21);                                                                 % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
fill(XB,YB,'k');                                                            % Plot polygon
if (strcmpi(load,'Circle'))                                                 % If circle is selected
    plot(x,y,'k--');                                                        % Plot dashed actual circle
end
for i = 1:1:numPan                                                          % Loop over all panels
    X(1) = XC(i);                                                           % Panel starting X point
    X(2) = XC(i) + S(i)*cosd(deltaD(i));                                    % Panel ending X point
    Y(1) = YC(i);                                                           % Panel starting Y point
    Y(2) = YC(i) + S(i)*sind(deltaD(i));                                    % Panel ending Y point
    if (i == 1)                                                             % For first panel
        p1 = plot(X,Y,'b-','LineWidth',2);                                  % Plot first panel normal vector
    elseif (i == 2)                                                         % For second panel
        p2 = plot(X,Y,'g-','LineWidth',2);                                  % Plot second panel normal vector
    else                                                                    % For every other panel
        plot(X,Y,'r-','LineWidth',2);                                       % Plot panel normal vector
    end
end
legend([p1,p2],{'Panel 1','Panel 2'});                                      % Add legend
xlabel('X Units');                                                          % Set X-label
ylabel('Y Units');                                                          % Set Y-label
axis equal;                                                                 % Set axes equal
zoom reset;                                                                 % Reset zoom


