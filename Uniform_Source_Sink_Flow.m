% ELEMENTARY FLOW: UNIFORM FLOW + SOURCE/SINK FLOW
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

Vinf   = 1;                                                                 % Freestream velocity [arb]
alpha  = 0;                                                                 % Freestream anlge of attack [deg]
lambda = 2;                                                                 % Source/sink strength
X0     = 0;                                                                 % Source/sink X coordinate
Y0     = 0;                                                                 % Source/sink Y coordainte

%% CALCULATIONS

% Create grid
numX    = 50;                                                               % Number of X points
numY    = 50;                                                               % Number of Y points
X       = linspace(-10,10,numX)';                                           % X-point array
Y       = linspace(-10,10,numY)';                                           % Y-pont array
[XX,YY] = meshgrid(X,Y);                                                    % Create the meshgrid

% Solve for velocities
Vx = zeros(numX,numY);                                                      % Initialize X velocity component
Vy = zeros(numX,numY);                                                      % Initialize Y velocity component
r  = zeros(numX,numY);                                                      % Initialize radius
for i = 1:1:numX                                                            % Loop over X points
    for j = 1:1:numY                                                        % Loop over Y points
        x       = XX(i,j);                                                  % X coordinate
        y       = YY(i,j);                                                  % Y coordinate
        dx      = x - X0;                                                   % X distance from source/sink
        dy      = y - Y0;                                                   % Y distance from source/sink
        r       = sqrt(dx^2 + dy^2);                                        % Distance from source/sink
        Vx(i,j) = Vinf*cosd(alpha) + ((lambda*dx)/(2*pi*r^2));              % Compute X velocity component
        Vy(i,j) = Vinf*sind(alpha) + ((lambda*dy)/(2*pi*r^2));              % Compute Y velocity component
    end
end

%% COMPUTE CIRCULATION

a    = 5;                                                                   % Horizontal axis half-length
b    = 5;                                                                   % Vertical axis half-length
x0   = 0;                                                                   % Ellipse center X coordinate
y0   = 0;                                                                   % Ellipse center Y coordinate
numT = 100;                                                                 % Number of points along ellipse
[Gamma,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);    % Call circulation calculation
fprintf('Circulation: %g\n',Gamma);                                         % Display circulation result

%% PLOTTING

% Streamline starting points
numSL  = 100;                                                               % Number of streamlines
xStart = -10.*ones(numSL,1);                                                % Streamline starting X coordinates
yStart = linspace(-10,10,numSL);                                            % Streamline starting Y coordinates

% Quiver plot
figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
quiver(X,Y,Vx,Vy,'r');                                                      % Plot velocity vectors
streamline(XX,YY,Vx,Vy,xStart,yStart)                                       % Plot streamlines
plot(xC,yC,'b-','LineWidth',2);                                             % Plot ellipse
xlim([-6 6]);                                                               % Set X-limits
ylim([-6 6]);                                                               % Set Y-limits
xlabel('X Axis');                                                           % Set X-label
ylabel('Y Axis');                                                           % Set Y-label
axis('equal');                                                              % Set axes equal
zoom reset;                                                                 % Reset zoom
