% SOURCE PANEL METHOD - SINGLE AIRFOIL
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 01/01/19
% Updated: 01/01/19 - Copied code from SP_Circle.m
%                   - Works as expected for airfoil geometry
%          01/11/19 - Added streamline plotting
%          02/09/20 - Added DAT airfoil loading option with XFOIL function
% Notes  : This code is not optimized, but is instead written in such a way
%          that it is easy to follow along with my YouTube video derivations
% 
% Functions Needed:
% - XFOIL.m
% - COMPUTE_IJ_SPM.m
% - STREAMLINE_SPM.m
% - COMPUTE_CIRCULATION.m
% 
% Programs Needed:
% - xfoil.exe
% 
% Folder Needed:
% - Airfoil_DAT_Selig: folder containing all Selig-format airfoils
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
% - [6]: UIUC Airfoil Database: Download All Files using Python
%           Link: https://www.youtube.com/watch?v=nILo18DlqAo
% - [7]: Python code for downloading Selig airfoil DAT files
%           Link: http://www.joshtheengineer.com/2019/01/30/uiuc-airfoil-database-file-download/

% clear;
clear;
clc;

%% KNOWNS

% Airfoil loading flags
flagAirfoil.XFoilCreate = 1;                                                % Create specified NACA airfoil in XFOIL
flagAirfoil.XFoilLoad   = 0;                                                % Load Selig-format airfoil from directory

% User-defined knowns
Vinf = 1;                                                                   % Freestream velocity []
AoA  = 0;                                                                   % Angle of attack [deg]
NACA = '2412';                                                              % NACA airfoil to load [####]

% Convert angle of attack to radians
alpha = AoA*(pi/180);                                                       % Angle of attack [rad]

% Plotting flags
flagPlot = [1;          % Airfoil with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison (XFOIL vs. SPM)
            1;          % Airfoil streamlines
            1];         % Pressure coefficient contour

%% XFOIL - CREATE/LOAD AIRFOIL

% PPAR menu options
PPAR.N  = '170';                                                            % "Number of panel nodes"
PPAR.P  = '4';                                                              % "Panel bunching parameter"
PPAR.T  = '1';                                                              % "TE/LE panel density ratios"
PPAR.R  = '1';                                                              % "Refined area/LE panel density ratio"
PPAR.XT = '1 1';                                                            % "Top side refined area x/c limits"
PPAR.XB = '1 1';                                                            % "Bottom side refined area x/c limits"

% Call XFOIL function to obtain the following:
% - Airfoil coordinates
% - Pressure coefficient along airfoil surface
% - Lift, drag, and moment coefficients
[xFoilResults,success] = XFOIL(NACA,PPAR,AoA,flagAirfoil);                  % Get XFOIL results for prescribed airfoil
if (success == 0)                                                           % If user canceled airfoil dialog box
    return;                                                                 % Exit the program
end

% Separate out results from XFOIL function results
afName  = xFoilResults.afName;                                              % Airfoil name
xFoilX  = xFoilResults.X;                                                   % X-coordinate for Cp result
xFoilY  = xFoilResults.Y;                                                   % Y-coordinate for Cp result
xFoilCP = xFoilResults.CP;                                                  % Pressure coefficient
XB      = xFoilResults.XB;                                                  % Boundary point X-coordinate
YB      = xFoilResults.YB;                                                  % Boundary point Y-coordinate
xFoilCL = xFoilResults.CL;                                                  % Lift coefficient
xFoilCD = xFoilResults.CD;                                                  % Drag coefficient
xFoilCM = xFoilResults.CM;                                                  % Moment coefficient

% Number of boundary points and panels
numPts = length(XB);                                                        % Number of boundary points
numPan = numPts - 1;                                                        % Number of panels (control points)

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(numPan,1);                                                     % Initialize edge value array
for i = 1:1:numPan                                                          % Loop over all panels
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              % Compute edge values
end
sumEdge = sum(edge);                                                        % Sum all edge values

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
sumLambda = sum(lambda.*S);                                                 % Check sum of source panel strengths
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

%% COMPUTE LIFT AND DRAG

% Compute normal and axial force coefficients
CN = -Cp.*S.*sin(beta);                                                     % Normal force coefficient []
CA = -Cp.*S.*cos(beta);                                                     % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               % Decompose axial and normal to lift coefficient []
CD = sum(CN.*sind(AoA)) + sum(CA.*cosd(AoA));                               % Decompose axial and normal to drag coefficient []
CM = sum(Cp.*(XC-0.25).*S.*cos(phi));                                       % Moment coefficient []

% Print the results to the Command Window
fprintf('======= RESULTS =======\n');
fprintf('Lift Coefficient (CL)\n');
fprintf('\tSPM  : %2.8f\n',CL);
fprintf('\tXFOIL: %2.8f\n',xFoilCL);
fprintf('Drag Coefficient (CD)\n');
fprintf('\tSPM  : %2.8f\n',CD);
fprintf('\tXFOIL: %2.8f\n',xFoilCD);
fprintf('Moment Coefficient (CM)\n');
fprintf('\tSPM  : %2.4f\n',CM);
fprintf('\tXFOIL: %2.4f\n',xFoilCM);

%% COMPUTE STREAMLINES

if (flagPlot(5) == 1 || flagPlot(6) == 1)
    % Grid parameters
    nGridX = 100;                                                           % X-grid for streamlines and contours
    nGridY = 100;                                                           % Y-grid for streamlines and contours
    xVals  = [-0.5; 1.5];                                                   % X-grid extents [min, max]
    yVals  = [-0.5; 0.5];                                                   % Y-grid extents [min, max]
    
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
    for m = 1:1:nGridX                                                      % Loop over X grid points
        for n = 1:1:nGridY                                                  % Loop over Y grid points
            XP = XX(m,n);                                                   % Current iteration's X grid point
            YP = YY(m,n);                                                   % Current iteration's Y grid point
            [Mx,My] = STREAMLINE_SPM(XP,YP,XB,YB,phi,S);                    % Compute Mx and My geometric integrals
            
            [in,on] = inpolygon(XP,YP,XB,YB);                               % See if points are in or on the airfoil
            if (in == 1 || on == 1)                                         % If the grid point is in or on the airfoil
                Vx(m,n) = 0;                                                % Set X-velocity equal to zero
                Vy(m,n) = 0;                                                % Set Y-velocity equal to zero
            else                                                            % If the grid point is outside the airfoil
                Vx(m,n) = Vinf*cosd(AoA) + sum(lambda.*Mx./(2*pi));         % Compute X-velocity
                Vy(m,n) = Vinf*sind(AoA) + sum(lambda.*My./(2*pi));         % Compute Y-velocity
            end
        end
    end
    
    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector []
    CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
end

%% CIRCULATION AND SOURCE STRENGTH CHECK

if (flagPlot(5) == 1 || flagPlot(6) == 1)                                   % If we are plotting 5 or 6
    % Compute circulation
    a    = 0.75;                                                            % Ellipse horizontal half-length
    b    = 0.25;                                                            % Ellipse vertical half-length
    x0   = 0.5;                                                             % Ellipse center X-coordinate
    y0   = 0;                                                               % Ellipse center Y-coordinate
    numT = 5000;                                                            % Number of points on ellipse
    [Circulation,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(a,b,x0,y0,numT,...    % Compute circulation around ellipse
                                                        Vx,Vy,XX,YY);
    
    % Print values to Command Window
    fprintf('Sum of L   : %g\n',sumLambda);                                 % Print sum of source strengths
    fprintf('Circulation: %g\n',Circulation);                               % Print circulation
end

%% PLOTTING

% FIGURE: Airfoil with panel normal vectors
if (flagPlot(1) == 1)
    figure(1);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    fill(XB,YB,'k');                                                        % Plot airfoil
    for i = 1:1:numPan                                                      % Loop over all panels
        X(1) = XC(i);                                                       % Set X start of panel orientation vector
        X(2) = XC(i) + S(i)*cosd(betaD(i)+AoA);                             % Set X end of panel orientation vector
        Y(1) = YC(i);                                                       % Set Y start of panel orientation vector
        Y(2) = YC(i) + S(i)*sind(betaD(i)+AoA);                             % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',2);                                       % Plot panel normal vector
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
	xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Geometry with the following indicated:
% - Boundary pts, control pts, first panel, second panel
if (flagPlot(2) == 1)
    figure(2);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XB,YB,'k-','LineWidth',3);                                         % Plot airfoil panels
    p1 = plot([XB(1) XB(2)],[YB(1) YB(2)],'g-','LineWidth',2);              % Plot first panel
    p2 = plot([XB(2) XB(3)],[YB(2) YB(3)],'m-','LineWidth',2);              % Plot second panel
    pB = plot(XB,YB,'ko','MarkerFaceColor','k');                            % Plot boundary points (black circles)
    pC = plot(XC,YC,'ko','MarkerFaceColor','r');                            % Plot control points (red circles)
    legend([pB,pC,p1,p2],...                                                % Show legend
           {'Boundary','Control','First Panel','Second Panel'});
	xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Cp vectors at airfoil control points
if (flagPlot(3) == 1)
    figure(3);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    Cps = abs(Cp*0.25);                                                     % Scale and make positive all Cp values
    for i = 1:1:length(Cps)                                                 % Loop over all panels
        X(1) = XC(i);                                                       % Control point X-coordinate
        X(2) = XC(i) + Cps(i)*cosd(betaD(i)+AoA);                           % Ending X-value based on Cp magnitude
        Y(1) = YC(i);                                                       % Control point Y-coordinate
        Y(2) = YC(i) + Cps(i)*sind(betaD(i)+AoA);                           % Ending Y-value based on Cp magnitude
        
        if (Cp(i) < 0)                                                      % If pressure coefficient is negative
            p{1} = plot(X,Y,'r-','LineWidth',2);                                   % Plot as a red line
        elseif (Cp(i) >= 0)                                                 % If pressure coefficient is zero or positive
            p{2} = plot(X,Y,'b-','LineWidth',2);                                   % Plot as a blue line
        end
    end
    fill(XB,YB,'k');                                                        % Plot the airfoil as black polygon
	legend([p{1},p{2}],{'Negative Cp','Positive Cp'});                      % Show legend
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Cp vectors at airfoil surface panels

if (flagPlot(4) == 1)
    figure(4);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    midIndX = floor(length(xFoilCP)/2);                                     % Airfoil middle index for XFOIL data
    midIndS = floor(length(Cp)/2);                                          % Airfoil middle index for SPM data
    pXu = plot(xFoilX(1:midIndX),xFoilCP(1:midIndX),'b-','LineWidth',2);    % Plot Cp for upper surface of airfoil from XFOIL
    pXl = plot(xFoilX(midIndX+1:end),xFoilCP(midIndX+1:end),'r-',...        % Plot Cp for lower surface of airfoil from XFOIL
                    'LineWidth',2);
    pSu = plot(XC(1:midIndS),Cp(1:midIndS),'ks','MarkerFaceColor','r');     % Plot Cp for upper surface of airfoil from SPM
    pSl = plot(XC(midIndS+1:end),Cp(midIndS+1:end),'ks',...                 % Plot Cp for lower surface of airfoil from SPM
                    'MarkerFaceColor','b');
    legend([pXu,pXl,pSu,pSl],...                                            % Show legend
           {'XFOIL Upper','XFOIL Lower','SPM Upper','SPM Lower'});
    xlabel('X Coordinate');                                                 % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 1]);                                                            % Set X-axis limits
    ylim('auto');                                                           % Set Y-axis limits to auto
    set(gca,'Ydir','reverse')                                               % Reverse direction of Y-axis
    title(['Airfoil: ' xFoilResults.afName ', CL/CL\_X = ' num2str(CL) '/' num2str(xFoilCL)]);  % Set title
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Airfoil streamlines
if (flagPlot(5) == 1)
    figure(5);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    for i = 1:1:length(Ysl)                                                 % Loop over all Y streamline starting points
        sl = streamline(XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % Plot the streamline
        set(sl,'LineWidth',2);                                              % Set streamline line width
    end
    fill(XB,YB,'k');                                                        % Plot airfoil as black polygon
	xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
    ylim(yVals);                                                            % Set Y-axis limits
	axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient contour
if (flagPlot(6) == 1)
    figure(6);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    contourf(XX,YY,CpXY,100,'EdgeColor','none');                            % Plot Cp contour
    fill(XB,YB,'k');                                                        % Plot airfoil as black polygon
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
    ylim(yVals);                                                            % Set Y-axis limits
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end
