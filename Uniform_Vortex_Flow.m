% ELEMENTARY FLOW: UNIFORM FLOW + VORTEX FLOW
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

Vinf  = 1;                                                                  % Freestream velocity [arb]
alpha = 0;                                                                  % Freestream angle of attack [deg]
gamma = 10;                                                                 % Vortex strength (+: CW, -: CCW)
X0    = 0;                                                                  % Vortex X coordinate
Y0    = 0;                                                                  % Vortex Y coordinate

%% CALCULATIONS

% Create grid
numX    = 50;                                                               % Number of X points
numY    = 50;                                                               % Number of Y points
X       = linspace(-10,10,numX)';                                           % X-point array
Y       = linspace(-10,10,numY)';                                           % Y-point array
[XX,YY] = meshgrid(X,Y);                                                    % Create the meshgrid

% Solve for velocities
Vx = zeros(numX,numY);                                                      % Initialize X velocity component
Vy = zeros(numX,numY);                                                      % Initialize Y velocity component
r  = zeros(numX,numY);                                                      % Initialize radius
for i = 1:1:numX                                                            % Loop over X points
    for j = 1:1:numY                                                        % Loop over Y points
        x       = XX(i,j);                                                  % X coordinate
        y       = YY(i,j);                                                  % Y coordinate
        dx      = x - X0;                                                   % X distance from vortex
        dy      = y - Y0;                                                   % Y distance from vortex
        r       = sqrt(dx^2 + dy^2);                                        % Distance from vortex
        Vx(i,j) = Vinf*cosd(alpha) + (gamma*dy)/(2*pi*r^2);                 % Compute X velocity component
        Vy(i,j) = Vinf*sind(alpha) + (-gamma*dx)/(2*pi*r^2);                % Compute Y velocity component
    end
end

%% COMPUTE CIRCULATION

% First circulation ellipse
a    = 5;                                                                       % Horizontal axis half-length
b    = 5;                                                                       % Vertical axis half-length
x0   = 0;                                                                       % Ellipse center X coordinate
y0   = 0;                                                                       % Ellipse center Y coordinate
numT = 1000;                                                                    % Number of points along ellipse 
[Gamma1,xC1,yC1,VxC1,VyC1] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);   % Call circulation calculation
fprintf('Circulation1: %g\n',Gamma1);                                           % Display circulation result

% Second circulation ellipse
a    = 1.5;                                                                     % Horizontal axis half-length
b    = 1.5;                                                                     % Vertical axis half-length
x0   = 0;                                                                       % Ellipse center X coordinate
y0   = 3;                                                                       % Ellipse center Y coordinate
numT = 1000;                                                                    % Number of points along ellipse
[Gamma2,xC2,yC2,VxC2,VyC2] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);   % Call circulation calculation
fprintf('Circulation2: %g\n',Gamma2);                                           % Display circulation result

%% PLOTTING

% Streamline starting points
numSL  = 20;                                                                % Number of streamlines
xStart = -10.*ones(numSL,1);                                                % Streamline starting X coordinates
yStart = linspace(-10,10,numSL);                                            % Streamline starting Y coordinates

% Quiver plot
figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
quiver(X,Y,Vx,Vy,'r');                                                      % Plot velocity vectors
h = streamline(XX,YY,Vx,Vy,xStart,yStart);                                  % Plot streamlines
set(h,'Color','k');                                                         % Set streamline color
p1 = plot(xC1,yC1,'b-','LineWidth',2);                                      % Plot first ellipse
p2 = plot(xC2,yC2,'m-','LineWidth',2);                                      % Plot second ellipse
xlim([-6 6]);                                                               % Set X-limits
ylim([-6 6]);                                                               % Set Y-limits
xlabel('X Axis');                                                           % Set X-label
ylabel('Y Axis');                                                           % Set Y-label
legend([p1,p2],{'\Gamma = 30.003','\Gamma = 0.005'},'Location','SE');       % Plot legend
axis('equal');                                                              % Set axes to equal sizes
zoom reset;                                                                 % Reset zoom
