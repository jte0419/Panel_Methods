% CIRCULATION
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 03/03/19
% Updated: 03/03/19 - Started code
%                   - Works as expected

clear;
clc;

%% KNOWNS

N = 5;                                                                      % Number of X and Y grid points

X  = linspace(0,N-1,N)';                                                    % X grid points
Y  = linspace(0,N-1,N)';                                                    % Y grid points
Vx = rand(N,N);                                                             % Random X velocities
Vy = rand(N,N);                                                             % Random Y velocities

[XX,YY] = meshgrid(X,Y);                                                    % Create the meshgrid

% Create ellipse
a    = 1.75;                                                                % Horizontal half-length
b    = 1.25;                                                                % Vertical half-length
x0   = 2;                                                                   % Center X-point
y0   = 2;                                                                   % Center Y-point
numT = 50;                                                                  % Number of ellipse points

[Gamma,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,Vx,Vy,XX,YY);    % Compute circulation
fprintf('Gamma: %g\n',Gamma);                                               % Print circulation

%% PLOTTING

figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Background color
set(gca,'FontSize',12);                                                     % Font size
quiver(XX,YY,Vx,Vy,'r-','AutoScale','off');                                 % Grid velocity vectors
plot(XX,YY,'ko','MarkerFaceColor','k');                                     % Grid points
quiver(xC,yC,VxC,VyC,'b-','AutoScale','off');                               % Ellipse velocity vectors
plot(xC,yC,'ko','MarkerFaceColor','b');                                     % Ellips points
xlim([-1 N]);                                                               % X-axis limits
ylim([-1 N]);                                                               % Y-axis limits
xlabel('X Axis');                                                           % X-axis label
ylabel('Y Axis');                                                           % Y-axis label

