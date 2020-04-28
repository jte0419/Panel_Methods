% VORTEX PANEL METHOD - SINGLE AIRFOIL
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 12/09/18
% Updated: 12/09/18 - Started code (based on earlier versions)
%                   - Works as expected
%          04/12/20 - Updated header info
% Notes  : This code is not optimized, but is instead written in such a way
%          that it is easy to follow along with my YouTube video derivations
% 
% Functions Needed:
% - XFOIL.m
% - COMPUTE_KL_VPM.m
% - STREAMLINE_VPM.m
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
% - [2]: Normal Geometric Integral VPM, K(ij)
%           Link: https://www.youtube.com/watch?v=5lmIv2CUpoc
% - [3]: Tangential Geometric Integral VPM, L(ij)
%           Link: https://www.youtube.com/watch?v=IxWJzwIG_gY
% - [4]: Streamline Geometric Integral VPM, Nx(pj) and Ny(pj)
%           Link: https://www.youtube.com/watch?v=TBwBnW87hso
% - [5]: Solving the System of Equations (VPM)
%           Link: https://www.youtube.com/watch?v=j3ETHFBiYOg
% - [6]: How To Compute Circulation
%           Link: https://www.youtube.com/watch?v=b8EnhiSjL3o
% - [7]: UIUC Airfoil Database: Download All Files using Python
%           Link: https://www.youtube.com/watch?v=nILo18DlqAo
% - [8]: Python code for downloading Selig airfoil DAT files
%           Link: http://www.joshtheengineer.com/2019/01/30/uiuc-airfoil-database-file-download/

clear;
clc;

%% KNOWNS

% Flag to specify creating or loading airfoil
flagAirfoil.XFoilCreate = 1;                                                % Create specified NACA airfoil in XFOIL
flagAirfoil.XFoilLoad   = 0;                                                % Load Selig-format airfoil from directory

% User-defined knowns
Vinf = 1;                                                                   % Freestream velocity []  (just leave this at 1)
AoA  = 0;                                                                   % Angle of attack [deg]
NACA = '0012';                                                              % NACA airfoil to load [####(#)]

% Convert angle of attack to radians
alpha = AoA*(pi/180);                                                       % Angle of attack [rad]

% Plotting flags
flagPlot = [0;          % Airfoil with panel normal vectors
            0;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison (XFOIL vs. VPM)
            0;          % Airfoil streamlines
            0];         % Pressure coefficient contour

%% XFOIL - CREATE/LOAD AIRFOIL

% PPAR menu options
PPAR.N  = '170';                                                            % "Number of panel nodes"
PPAR.P  = '4';                                                              % "Panel bunching parameter"
PPAR.T  = '1.5';                                                            % "TE/LE panel density ratios"
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

% Find geometric quantities of the airfoil
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face) [deg]
    if (phiD(i) < 0)                                                        % Make all panel angles positive [deg]
        phiD(i) = phiD(i) + 360;
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make all panel angles between 0 and 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

%% COMPUTE VORTEX PANEL STRENGTHS - REF [5]

% Geometric integral (normal [K] and tangential [L])
% - Refs [2] and [3]
[K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S);                                  % Compute geometric integrals

% Populate A matrix
% - Simpler option: A = -K;
A = zeros(numPan,numPan);                                                   % Initialize the A matrix
for i = 1:1:numPan                                                          % Loop over all i panels
    for j = 1:1:numPan                                                      % Loop over all j panels
        if (j == i)                                                         % If the panels are the same
            A(i,j) = 0;                                                     % Set A equal to zero
        else                                                                % If panels are not the same
            A(i,j) = -K(i,j);                                               % Set A equal to negative geometric integral
        end
    end
end

% Populate b array
% - Simpler option: b = -Vinf*2*pi*cos(beta);
b = zeros(numPan,1);                                                        % Initialize the b array
for i = 1:1:numPan                                                          % Loop over all panels
    b(i) = -Vinf*2*pi*cos(beta(i));                                         % Compute RHS array
end

% Satisfy the Kutta condition
pct    = 100;                                                               % Panel replacement percentage
panRep = floor((pct/100)*numPan);                                           % Replace this panel with Kutta condition eqn.
if (panRep == 0)                                                            % If we specify the first panel
    panRep = 1;                                                             % Make sure the index is not zero
end
A(panRep,:)   = 0;                                                          % Set all columns of the replaced panel equal to zero
A(panRep,1)   = 1;                                                          % Set first column of replaced panel equal to 1
A(panRep,end) = 1;                                                          % Set last column of replaced panel equal to 1
b(panRep)     = 0;                                                          % Set replaced panel value in b array equal to zero

% Compute gamma values from system of equations
gamma = A\b;                                                                % Compute all vortex strength values

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities on each panel
Vt = zeros(numPan,1);                                                       % Initialize tangential velocity array
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient array
for i = 1:1:numPan                                                          % Loop over all i panels
    addVal  = 0;                                                            % Reset the summation value to zero
    for j = 1:1:numPan                                                      % Loop over all j panels
        addVal = addVal - (gamma(j)/(2*pi))*L(i,j);                         % Sum all tangential vortex panel terms
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal + gamma(i)/2;                        % Compute tangential velocity by adding uniform flow and i=j terms
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient
end

%% COMPUTE LIFT AND MOMENT COEFFICIENTS

% Compute normal and axial force coefficients
CN = -Cp.*S.*sin(beta);                                                     % Normal force coefficient []
CA = -Cp.*S.*cos(beta);                                                     % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               % Decompose axial and normal to lift coefficient []
CM = sum(Cp.*(XC-0.25).*S.*cos(phi));                                       % Moment coefficient []

% Print the results to the Command Window
fprintf('======= RESULTS =======\n');
fprintf('Lift Coefficient (CL)\n');
fprintf('\tK-J  : %g\n',2*sum(gamma.*S));                                   % From Kutta-Joukowski lift equation
fprintf('\tVPM  : %2.8f\n',CL);                                             % From this VPM code
fprintf('\tXFOIL: %2.8f\n',xFoilCL);                                        % From XFOIL program
fprintf('Moment Coefficient (CM)\n');
fprintf('\tVPM  : %2.4f\n',CM);                                             % From this VPM code
fprintf('\tXFOIL: %2.4f\n',xFoilCM);                                        % From XFOIL program

%% COMPUTE STREAMLINES - REF [4]

if (flagPlot(5) == 1 || flagPlot(6) == 1)                                   % If we are plotting 5 or 6
    % Grid parameters
    nGridX = 150;                                                           % X-grid for streamlines and contours
    nGridY = 150;                                                           % Y-grid for streamlines and contours
	xVals  = [-0.5; 1.5];                                                   % X-grid extents [min, max]
    yVals  = [-0.3; 0.3];                                                   % Y-grid extents [min, max]
    
    % Streamline parameters
    stepsize = 0.01;                                                        % Step size for streamline propagation
    maxVert  = nGridX*nGridY*100;                                           % Maximum vertices
    slPct    = 25;                                                          % Percentage of streamlines of the grid
    Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points
    
    % Generate the grid points
    Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid
    Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid
    [XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays
    
    % Initialize velocities
    Vx = zeros(nGridX,nGridY);                                              % Initialize X-velocity matrix
    Vy = zeros(nGridX,nGridY);                                              % Initialize Y-velocity matrix
    
    % Solve for grid point X and Y velocities
    for m = 1:1:nGridX
        for n = 1:1:nGridY
            XP = XX(m,n);                                                   % Current iteration's X grid point
            YP = YY(m,n);                                                   % Current iteration's Y grid point
            [Nx,Ny] = STREAMLINE_VPM(XP,YP,XB,YB,phi,S);                    % Compute Nx and Ny geometric integrals
            
            [in,on] = inpolygon(XP,YP,XB,YB);
            if (in == 1 || on == 1)                                         % If the grid point is in or on the airfoil
                Vx(m,n) = 0;                                                % Set X-velocity equal to zero
                Vy(m,n) = 0;                                                % Set Y-velocity equal to zero
            else                                                            % If the grid point is outside the airfoil
                Vx(m,n) = Vinf*cosd(AoA) + sum(-gamma.*Nx./(2*pi));         % Compute X-velocity
                Vy(m,n) = Vinf*sind(AoA) + sum(-gamma.*Ny./(2*pi));         % Compute Y-velocity
            end
        end
    end
    
    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector []
    CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
end

%% CIRCULATION AND VORTEX STRENGTH CHECK - REF [6]

if (flagPlot(5) == 1 || flagPlot(6) == 1)                                   % If we are plotting streamlines or Cp contours
    % Compute circulation
    aa   = 0.75;                                                            % Ellipse horizontal half-length
    bb   = 0.25;                                                            % Ellipse vertical half-length
    x0   = 0.5;                                                             % Ellipse center X-coordinate
    y0   = 0;                                                               % Ellipse center Y-coordinate
    numT = 5000;                                                            % Number of points on ellipse
    [Circulation,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(aa,bb,x0,y0,numT,...  % Compute circulation around ellipse
                                                        Vx,Vy,XX,YY);
    
    % Print values to Command Window
    fprintf('======= CIRCULATION RESULTS =======\n');
    fprintf('Sum of G   : %g\n',sum(gamma.*S));                             % Print sum of vortex strengths
    fprintf('Circulation: %g\n',Circulation);                               % Print circulation
    fprintf('K-J Lift   : %g\n',2*Circulation);                             % Lift coefficient from K-J equation
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

% FIGURE: Pressure coefficient comparison (XFOIL vs. VPM)
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
    pVl = plot(XC(1:midIndS),Cp(1:midIndS),'ks','MarkerFaceColor','r');     % Plot Cp for upper surface of airfoil from SPM
    pVu = plot(XC(midIndS+1:end),Cp(midIndS+1:end),'ks',...                 % Plot Cp for lower surface of airfoil from SPM
                    'MarkerFaceColor','b');
    legend([pXu,pXl,pVu,pVl],...                                            % Show legend
           {'XFOIL Upper','XFOIL Lower','VPM Upper','VPM Lower'});
    xlabel('X Coordinate');                                                 % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 1]);                                                            % Set X-axis limits
    ylim('auto');                                                           % Set Y-axis limits to auto
    set(gca,'Ydir','reverse')                                               % Reverse direction of Y-axis
    title(['Airfoil: ' xFoilResults.afName ...                              % Title
           ', CL_{VPM}/CL_{XFOIL} = ' ...
           num2str(CL,4) '/' num2str(xFoilCL,4)]);
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
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
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
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
end
