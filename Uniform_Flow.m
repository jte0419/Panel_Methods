% ELEMENTARY FLOW - UNIFORM FLOW
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
alpha = 0;                                                                  % Angle of attack [deg]

%% CALCULATIONS

% Create grid
numX    = 10;                                                               % Number of X points
numY    = 10;                                                               % Number of Y points
X       = linspace(-10,10,numX)';                                           % X-point array
Y       = linspace(-10,10,numY)';                                           % Y-point array
[XX,YY] = meshgrid(X,Y);                                                    % Create the meshgrid

% Solve for velocities
Vx = zeros(numX,numY);                                                      % Initialize X velocity component
Vy = zeros(numX,numY);                                                      % Initialize Y velocity component
for i = 1:1:numX                                                            % Loop over X points
    for j = 1:1:numY                                                        % Loop over Y points
        Vx(i,j) = Vinf*cosd(alpha);                                         % Compute X velocity component
        Vy(i,j) = Vinf*sind(alpha);                                         % Compute Y velocity component
    end
end

%% COMPUTE CIRCULATION

a    = 5;                                                                   % Horizontal axis half-length
b    = 5;                                                                   % Vertical axis half-length
x0   = 0;                                                                   % Ellipse center X coordinate
y0   = 0;                                                                   % Ellipse center Y coordinate
numT = 100;                                                                 % Number of points along ellipse
[Gamma,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);    % Call circulation calculation function
fprintf('Circulation: %g\n',Gamma);                                         % Display circulation result

%% PLOTTING

% Quiver plot
figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
quiver(X,Y,Vx,Vy,'r');                                                      % Plot velocity vectors
plot(xC,yC,'b-','LineWidth',2);                                             % Plot ellipse
xlim([min(X) max(X)]);                                                      % Set X-limits
ylim([min(Y) max(Y)]);                                                      % Set Y-limits
xlabel('X Axis');                                                           % Set X-label
ylabel('Y Axis');                                                           % Set Y-label
title('Uniform Flow');                                                      % Set title
axis('equal');                                                              % Set axes to equal sizes
zoom reset;                                                                 % Reset zoom
