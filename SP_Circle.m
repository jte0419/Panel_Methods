% SOURCE PANEL METHOD - CIRCLE
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 01/01/19
% Updated: 01/01/19 - Started code
%                   - Works as expected
%          01/11/19 - Added streamline plotting
% Notes  : This code is not optimized, but is instead written in such a way
%          that it is easy to follow along with my YouTube video derivations
% 
% Functions Needed:
% - COMPUTE_IJ_SPM.m
% - STREAMLINE_SPM.m
% 
% References
% - [1]: Panel Method Geometry
%           Link: https://www.youtube.com/watch?v=kIqxbd937PI
% - [2]: Normal Geometric Integral SPM, I(ij)
%           Link: https://www.youtube.com/watch?v=76vPudNET6U
% - [3]: Tangential Geometric Integral SPM, J(ij)
%           Link: https://www.youtube.com/watch?v=JRHnOsueic8
% - [4]: Streamline Geometric Integral SPM, Mx(ij) and My(ij)
%           Link: https://www.youtube.com/watch?v=BnPZjGCatcg
% - [5]: Solving the System of Equations
%           Link: https://www.youtube.com/watch?v=ep7vPzGYsbw

clear;
clc;

%% KNOWNS

% User-defined knowns
Vinf = 1;                                                                   % Freestream velocity
AoA  = 0;                                                                   % Angle of attack [deg]
numB = 9;                                                                   % Number of boundary points (including endpoint)
tO   = (360/(numB-1))/2;                                                    % Boundary point angle offset [deg]

% Plotting flags
flagPlot = [1;          % Shape polygon with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Analytical and SPM pressure coefficient plot
            1;          % Streamlines
            1];         % Pressure coefficient contours

%% CREATE CIRCLE BOUNDARY POINTS

% Angles used to compute boundary points
theta = linspace(0,360,numB)';                                              % Create angles for computing boundary point locations [deg]
theta = theta + tO;                                                         % Add panel angle offset [deg]
theta = theta*(pi/180);                                                     % Convert from degrees to radians [rad]

% Boundary points
XB = cos(theta);                                                            % Compute boundary point X-coordinate (radius of 1)
YB = sin(theta);                                                            % Compute boundary point Y-coordinate (radius of 1)

% Number of panels
numPan = length(XB)-1;                                                      % Number of panels (control points)

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(numPan,1);                                                     % Initialize edge value array
for i = 1:1:numPan                                                          % Loop over all panels
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              % Compute edge value
end
sumEdge = sum(edge);                                                        % Sum of all edge values

% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            % If panels are CCW
    XB = flipud(XB);                                                        % Flip the X-data array
    YB = flipud(YB);                                                        % Flip the Y-data array
end

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
XC   = zeros(numPan,1);                                                     % Initialize control point X-coordinate array
YC   = zeros(numPan,1);                                                     % Initialize control point Y-coordinate array
S    = zeros(numPan,1);                                                     % Initialize panel length array
phiD = zeros(numPan,1);                                                     % Initialize panel orientation angle array [deg]

% Find geometric quantities of airfoil
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face)
    if (phiD(i) < 0)                                                        % Make all panel angles positive
        phiD(i) = phiD(i) + 360;
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make sure angles aren't greater than 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

%% COMPUTE SOURCE PANEL STRENGTHS - REF [5]

% Geometric integral (normal [I] and tangential [J])
% - Refs [2] and [3]
[I,J] = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S);                                  % Compute geometric integrals

% Populate A matrix
% - Simpler option: A = I + pi*eye(numPan,numPan);
A = zeros(numPan,numPan);                                                   % Initialize the A matrix
for i = 1:1:numPan                                                          % Loop over all i panels
    for j = 1:1:numPan                                                      % Loop over all j panels
        if (i == j)                                                         % If the panels are the same
            A(i,j) = pi;                                                    % Set A equal to pi
        else                                                                % If panels are not the same
            A(i,j) = I(i,j);                                                % Set A equal to geometric integral
        end
    end
end

% Populate b array
% - Simpler option: b = -Vinf*2*pi*cos(beta);
b = zeros(numPan,1);                                                        % Initialize the b array
for i = 1:1:numPan                                                          % Loop over all panels
    b(i) = -Vinf*2*pi*cos(beta(i));                                         % Compute RHS array
end

% Compute source panel strengths (lambda) from system of equations
lambda  = A\b;                                                              % Compute all source strength values

% Check the sum of the source strenghts
% - This should be very close to zero for a closed polygon
fprintf('Sum of L: %g\n',sum(lambda.*S));                                   % Print sum of all source strengths

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities
% - Simpler method: Vt = Vinf*sin(beta) + J*lambda/(2*pi);
%                   Cp = 1-(Vt/Vinf).^2;
Vt = zeros(numPan,1);                                                       % Initialize tangential velocity array
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient array
for i = 1:1:numPan                                                          % Loop over all i panels
    addVal  = 0;                                                            % Reset the summation value to zero
    for j = 1:1:numPan                                                      % Loop over all j panels
        addVal = addVal + (lambda(j)/(2*pi))*(J(i,j));                      % Sum all tangential source panel terms
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal;                                     % Compute tangential velocity by adding uniform flow term
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient
end

% Analytical angles and pressure coefficients
analyticTheta = linspace(0,2*pi,200)';                                      % Analytical theta angles [rad]
analyticCP    = 1-4*sin(analyticTheta).^2;                                  % Analytical pressure coefficient []

%% COMPUTE LIFT AND DRAG

% Compute normal and axial force coefficients
CN = -Cp.*S.*sin(beta);                                                      % Normal force coefficient []
CA = -Cp.*S.*cos(beta);                                                      % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               % Decompose axial and normal to lift coefficient []
CD = sum(CN.*sind(AoA)) + sum(CA.*cosd(AoA));                               % Decompose axial and normal to drag coefficient []

% Display lift and drag coefficients in command window
fprintf('CL      : %g\n',CL);                                               % Display lift coefficient (should be zero)
fprintf('CD      : %g\n',CD);                                               % Display drag coefficient (should be zero)

%% COMPUTE STREAMLINES - REF [4]

if (flagPlot(4) == 1 || flagPlot(5) == 1)
    % Grid parameters
    nGridX = 100;                                                           % X-grid for streamlines and contours
    nGridY = 100;                                                           % Y-grid for streamlines and contours
    xVals  = [-1.5; 1.5];                                                   % X-grid extents [min, max]
    yVals  = [-1.5; 1.5];                                                   % Y-grid extents [min, max]

    % Streamline parameters
    stepsize = 0.01;                                                        % Step size for streamline propagation
    maxVert  = nGridX*nGridY*10;                                            % Maximum vertices
    slPct    = 30;                                                          % Percentage of streamlines of the grid
    Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points

    % Generate the grid points
    Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid
    Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid
    [XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays

    % Initialize velocities
    Vx = zeros(nGridX,nGridY);                                              % Initialize X velocity matrix
    Vy = zeros(nGridX,nGridY);                                              % Initialize Y velocity matrix

    % Solve for grid point X and Y velocities
    for m = 1:1:nGridX                                                      % Loop over X-grid points
        for n = 1:1:nGridY                                                  % Loop over Y-grid points
            XP = XX(m,n);                                                   % Isolate X point
            YP = YY(m,n);                                                   % Isolate Y point
            [Mx,My] = STREAMLINE_SPM(XP,YP,XB,YB,phi,S);                    % Compute streamline Mx and My values (Ref [4])

            % Check if grid points are in object
            % - If they are, assign a velocity of zero
            [in,on] = inpolygon(XP,YP,XB,YB);                               % Find whether (XP,YP) is in polygon body
            if (in == 1 || on == 1)                                         % If (XP, YP) is in the polygon body
                Vx(m,n) = 0;                                                % X-velocity is zero
                Vy(m,n) = 0;                                                % Y-velocity is zero
            else                                                            % If (XP,YP) is not in the polygon body
                Vx(m,n) = Vinf*cosd(AoA) + sum(lambda.*Mx./(2*pi));         % Compute X-velocity
                Vy(m,n) = Vinf*sind(AoA) + sum(lambda.*My./(2*pi));         % Compute Y-velocity
            end
        end
    end

    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector
    CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
end

%% PLOTTING

% FIGURE: Shape polygon with panel normal vectors
if (flagPlot(1) == 1)
    figure(1);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    fill(XB,YB,'k');                                                        % Plot polygon
    for i = 1:1:numPan                                                      % Loop over all panels
        X(1) = XC(i);                                                       % Set X start of panel orientation vector
        X(2) = XC(i) + S(i)*cosd(betaD(i)+AoA);                             % Set X end of panel orientation vector
        Y(1) = YC(i);                                                       % Set Y start of panel orientation vector
        Y(2) = YC(i) + S(i)*sind(betaD(i)+AoA);                             % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',3);                                       % Plot panel normal vector
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Geometry with the following indicated:
% - Boundary points, control points, first panel, second panel
if (flagPlot(2) == 1)
    figure(2);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XB,YB,'k-','LineWidth',3);                                         % Plot panels
    p1 = plot([XB(1) XB(2)],[YB(1) YB(2)],'g-','LineWidth',3);              % Plot first panel
    p2 = plot([XB(2) XB(3)],[YB(2) YB(3)],'m-','LineWidth',3);              % Plot second panel
    pB = plot(XB,YB,'ko','MarkerFaceColor','k','MarkerSize',10);            % Plot boundary points
    pC = plot(XC,YC,'ko','MarkerFaceColor','r','MarkerSize',10);            % Plot control points
    legend([pB,pC,p1,p2],...                                                % Show legend
           {'Boundary','Control','First Panel','Second Panel'});
    ylim([-1 1]);                                                           % Set Y-limits
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Analytical and SPM pressure coefficient plot
if (flagPlot(3) == 1)
    figure(3);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    pA = plot(analyticTheta,analyticCP,'k-','LineWidth',3);                 % Plot analytical pressure coefficient
    pC = plot(beta,Cp,'ks','MarkerFaceColor','r','MarkerSize',10);          % Plot compute pressure coefficient
    xlabel('Angle [rad]');                                                  % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 2*pi]);                                                         % Set X-limits
    ylim([-3.5 1.5]);                                                       % Set Y-limits
    legend([pA,pC],{'Analytical','SPM'},'Location','S');                    % Add legend
end

% FIGURE: Streamlines (and quiver if commented in)
if (flagPlot(4) == 1)
    figure(4);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
%     quiver(XX,YY,Vx,Vy,'r');                                                % Plot velocity vectors
    for i = 1:1:length(Ysl)                                                 % Loop over all streamlines
        sl = streamline(XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % streamline(X,Y,U,V,startx,starty)
        set(sl,'LineWidth',3);                                              % Change streamline width
    end
    fill(XB,YB,'k');                                                        % Plot polygon
    xlim(xVals);                                                            % Set X-limits
    ylim(yVals);                                                            % Set Y-limits
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient contours
if (flagPlot(5) == 1)
    figure(5);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    contourf(XX,YY,CpXY,100,'EdgeColor','none');                            % Plot contour
    fill(XB,YB,'k');                                                        % Plot polygon
    xlim(xVals);                                                            % Set X-limits
    ylim(yVals);                                                            % Set Y-limits
    zoom reset;                                                             % Reset zoom
end
