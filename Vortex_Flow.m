% ELEMENTARY FLOW: VORTEX FLOW
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 02/09/19
% Updated: 02/09/19 - Started code
%                   - Works as expected
%          02/11/19 - Added circulation calculation

clear;
clc;

%% KNOWNS

gamma = 1;                                                                  % Vortex strength (+: CW, -: CCW)
X0    = 0;                                                                  % Vortex X coordinate
Y0    = 0;                                                                  % Vortex Y coordinate

%% CALCULATIONS

% Create grid
numX    = 100;                                                              % Number of X points
numY    = 100;                                                              % Number of Y points
X       = linspace(-10,10,numX)';                                           % X-point array
Y       = linspace(-10,10,numY)';                                           % Y-point array
[XX,YY] = meshgrid(X,Y);                                                    % Create the meshgrid

% Solve for velocities
Vx = zeros(numX,numY);                                                      % Initialize X velocity component
Vy = zeros(numX,numY);                                                      % Initialize Y velocity component
V  = zeros(numX,numY);                                                      % Initialize velocity
Vt = zeros(numX,numY);                                                      % Initialize tangential velocity component
Vr = zeros(numX,numY);                                                      % Initialize radial velocity component
r  = zeros(numX,numY);                                                      % Intialize radius
for i = 1:1:numX                                                            % Loop over X points
    for j = 1:1:numY                                                        % Loop over Y points
        x       = XX(i,j);                                                  % X coordinate
        y       = YY(i,j);                                                  % Y coordinate
        dx      = x - X0;                                                   % X distance from vortex
        dy      = y - Y0;                                                   % Y distance from vortex
        r       = sqrt(dx^2 + dy^2);                                        % Distance from vortex
        Vx(i,j) = (gamma*dy)/(2*pi*r^2);                                    % Compute X velocity component
        Vy(i,j) = (-gamma*dx)/(2*pi*r^2);                                   % Compute Y velocity component
        V(i,j) = sqrt(Vx(i,j)^2 + Vy(i,j)^2);                               % Compute velocity
        Vt(i,j) = -gamma/(2*pi*r);                                          % Compute tangential velocity component
        Vr(i,j) = 0;                                                        % Compute radial velocity component
    end
end

%% COMPUTE CIRCULATION

a    = 2;                                                                   % Horizontal axis half-length
b    = 2;                                                                   % Vertical axis half-length
x0   = 0;                                                                   % Ellipse center X coordinate
y0   = 0;                                                                   % Ellipse center Y coordinate
numT = 100;                                                                 % Number of points along ellipse
[Gamma,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);    % Call circulation calculation
fprintf('Circulation: %g\n',Gamma);                                         % Display circulation result

%% PLOTTING

% Plot velocity vectors
figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
quiver(X,Y,Vx,Vy,'r');                                                      % Plot velocity vectors
plot(xC,yC,'b-','LineWidth',2);                                             % Plot ellipse
xlim([-3 3]);                                                               % Set X-limits
ylim([-3 3]);                                                               % Set Y-limits
xlabel('X Axis');                                                           % Set X-label
ylabel('Y Axis');                                                           % Set Y-label
axis('equal');                                                              % Set axes to equal sizes
zoom reset;                                                                 % Reset zoom
