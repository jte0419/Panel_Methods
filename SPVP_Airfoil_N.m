% SOURCE/VORTEX PANEL METHOD - MULTI-AIRFOIL
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started: 11/14/19
% Updated: 11/14/19 - Started code (based on earlier versions)
%          11/15/19 - Works as expected
%          04/12/20 - Updated header info
%          04/16/20 - Adding multi-element code from other file (01/21/19)
%          04/18/20 - Updating SPVP code for multi-element
% Note 1 : This code is not optimized, but is instead written in such a way
%          that it is easy to follow along with my YouTube video derivations
% Note 2 : Don't use this code for anything serious
% 
% Functions Needed:
% - XFOIL_N.m
% - COMPUTE_IJ_SPM_N.m
% - COMPUTE_KL_VPM_N.m
% - STREAMLINE_SPM_N.m
% - STREAMLINE_VPM_N.m
% - COMPUTE_CIRCULATION.m
% 
% Programs Needed:
% - xfoil.exe
% 
% Folder Needed:
% - Airfoil_DAT_Selig: folder containing all Selig-format airfoils
% 
% References
% - [1] : Panel Method Geometry
%           Link: https://www.youtube.com/watch?v=kIqxbd937PI
% - [2] : Normal Geometric Integral SPM, I(ij)
%           Link: https://www.youtube.com/watch?v=76vPudNET6U
% - [3] : Tangential Geometric Integral SPM, J(ij)
%           Link: https://www.youtube.com/watch?v=JRHnOsueic8
% - [4] : Streamline Geometric Integral SPM, Mx(pj) and My(pj)
%           Link: https://www.youtube.com/watch?v=BnPZjGCatcg
% - [5] : Solving the System of Equations (SPM)
%           Link: https://www.youtube.com/watch?v=ep7vPzGYsbw
% - [6] : Normal Geometric Integral VPM, K(ij)
%           Link: https://www.youtube.com/watch?v=5lmIv2CUpoc
% - [7] : Tangential Geometric Integral VPM, L(ij)
%           Link: https://www.youtube.com/watch?v=IxWJzwIG_gY
% - [8] : Streamline Geometric Integral VPM, Nx(pj) and Ny(pj)
%           Link: https://www.youtube.com/watch?v=TBwBnW87hso
% - [9] : Solving the System of Equations (VPM)
%           Link: https://www.youtube.com/watch?v=j3ETHFBiYOg
% - [10]: Source/Vortex Panel Method System of Equations
%           Link: https://www.youtube.com/watch?v=bc_pkKGEypU
% - [11]: Source/Vortex Panel Method: Airfoil
%           Link: https://www.youtube.com/watch?v=V77QTAgZuqw
% - [12]: How To Compute Circulation
%           Link: https://www.youtube.com/watch?v=b8EnhiSjL3o
% - [13]: UIUC Airfoil Database: Download All Files using Python
%           Link: https://www.youtube.com/watch?v=nILo18DlqAo
% - [14]: Python code for downloading Selig airfoil DAT files
%           Link: http://www.joshtheengineer.com/2019/01/30/uiuc-airfoil-database-file-download/
% 
% Some New Variable Definitions:
% - numAF     : Total number of airfoils defined by the user [#]
% - numPtsAF  : Number of boundary points for each airfoil [numAF x 1]
% - numPanAF  : Number of panels for each airfoil [numAF x 1]
% - totPan    : Number of all panels, including false panels [#]
% - falsePan  : Indices of the false panels [numAF-1 x 1]
% - panInd    : Indices of the real panels [length(numPanAF) x 1]
%               iInd and jInd are the same as panInd, but for use in functions
% - kS/kE     : Starting/ending indices of the airfoils [numAF x 1]
% - numPtsTot : Total number of boundary points [#]
% - numPanTot : Total number of panels [#]
% - numPan    : Total number of real airfoil panels [#]

clear;
clc;

%% KNOWNS

% % Triple Airfoil
% % - Used for explanation of code updates
% AF_name   = {'0012'; '0012'; '0012'};                                       % NACA airfoils to load [####(#)]
% AF_load   = [0; 0; 0];                                                      % Load [1] or create [0] airfoil
% numPtsTot = {'5'; '5'; '5'};                                                % Number of points per airfoil
% AF_flip   = [1; 1; 1];                                                      % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 1; 1];                                                      % Scale airfoil
% AF_angle  = [0; 5; 20];                                                     % Angle of rotation of airfoil [deg]
% AF_offset = [0     0;                                                       % Offset from origin [x y]
%              1.5   -0.3;
%              3.0   -0.6];

% % Single Airfoil
% % - Tests symmetric single airfoil
% AF_name   = {'0012'};                                                       % NACA airfoils to load [####(#)]
% AF_load   = [0];                                                            % Load [1] or create [0] airfoil
% numPtsTot = {'170'};                                                        % Number of points per airfoil
% AF_flip   = [1];                                                            % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1];                                                            % Scale airfoil
% AF_angle  = [0];                                                            % Angle of rotation of airfoil [deg]
% AF_offset = [0 0];                                                          % Offset from origin [x y]

% % Single Airfoil
% % - Tests single cambered airfoil, compare to SPVP_Airfoil.m
% AF_name   = {'2412'};                                                       % NACA airfoils to load [####(#)]
% AF_load   = [0];                                                            % Load [1] or create [0] airfoil
% numPtsTot = {'170'};                                                        % Number of points per airfoil
% AF_flip   = [1];                                                            % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1];                                                            % Scale airfoil
% AF_angle  = [0];                                                            % Angle of rotation of airfoil [deg]
% AF_offset = [0 0];                                                          % Offset from origin [x y]

% % Double Airfoil
% % - Tests double symmetric airfoil at 0 deg AoA
% AF_name   = {'0012'; '0012'};                                               % NACA airfoils to load [####(#)]
% AF_load   = [0; 0];                                                         % Load [1] or create [0] airfoil
% numPtsTot = {'100'; '100'};                                                 % Number of points per airfoil
% AF_flip   = [1; 1];                                                         % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 1];                                                         % Scale airfoil
% AF_angle  = [0; 0];                                                         % Angle of rotation of airfoil [deg]
% AF_offset = [0   0;                                                         % Offset from origin [x y]
%              1.2 0];

% % Triple Airfoil
% % - Tests triple symmetric airfoil at 0 deg AoA
% AF_name   = {'0012'; '0012'; '0012'};                                       % NACA airfoils to load [####(#)]
% AF_load   = [0; 0; 0];                                                      % Load [1] or create [0] airfoil
% numPtsTot = {'35'; '35'; '35'};                                             % Number of points per airfoil
% AF_flip   = [1; 1; 1];                                                      % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 1; 1];                                                      % Scale compared to first airfoil
% AF_angle  = [0; 0; 0];                                                      % Angle of rotation of airfoil [deg]
% AF_offset = [0     0;                                                       % Offset from origin [x y]
%              1.2   0;
%              2.4   0];

% % Double Airfoil
% % - Tests loading and creating of airfoils in same airfoil system
% AF_name   = {'nlr7301'; '2412'};                                            % NACA airfoils to load [####(#)]
% AF_load   = [1; 0];                                                         % Load [1] or create [0] airfoil
% numPtsTot = {'100'; '100'};                                                 % Number of points per airfoil
% AF_flip   = [1; 1];                                                         % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 1];                                                         % Scale airfoil
% AF_angle  = [0; 0];                                                         % Angle of rotation of airfoil [deg]
% AF_offset = [0   0;                                                         % Offset from origin [x y]
%              1.2 0];

% % Double Airfoil (F1 Rear Wing, DRS Closed)
% % - Tests inversion for racecar rear wing (no ground effect)
% AF_name   = {'8412'; '5412'};                                               % NACA airfoils to load [####(#)]
% AF_load   = [0; 0];                                                         % Load [1] or create [0] airfoil
% numPtsTot = {'100'; '100'};                                                 % Number of points per airfoil
% AF_flip   = [-1; -1];                                                       % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 0.35];                                                      % Scale airfoil
% AF_angle  = [-15; -40];                                                     % Angle of rotation of airfoil [deg]
% AF_offset = [0   0;                                                         % Offset from origin [x y]
%              0.9 0.3];

% % Double Airfoil (F1 Rear Wing, DRS Open)
% % - Tests inversion for racecar rear wing (no ground effect)
% AF_name   = {'8412'; '5412'};                                               % NACA airfoils to load [####(#)]
% AF_load   = [0; 0];                                                         % Load [1] or create [0] airfoil
% numPtsTot = {'100'; '100'};                                                 % Number of points per airfoil
% AF_flip   = [-1; -1];                                                       % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 0.35];                                                      % Scale airfoil
% AF_angle  = [-15; -15];                                                     % Angle of rotation of airfoil [deg]
% AF_offset = [0   0;                                                         % Offset from origin [x y]
%              0.85 0.43];

% % Triple Airfoil
% % - Just a fun test case with three cambered airfoils
% AF_name   = {'6412'; '3412'; '2412'};                                       % NACA airfoils to load [####(#)]
% AF_load   = [0; 0; 0];                                                      % Load [1] or create [0] airfoil
% numPtsTot = {'35'; '35'; '35'};                                             % Number of points per airfoil
% AF_flip   = [1; 1; 1];                                                      % Flip the airfoil vertically [1 = NO, -1 = YES]
% AF_scale  = [1; 0.3; 0.2];                                                  % Scale airfoil
% AF_angle  = [0; 25;  45];                                                   % Angle of rotation of airfoil [deg]
% AF_offset = [0     0;                                                       % Offset from origin [x y]
%              1.05 -0.05;
%              1.35  -0.25];

% Narsipur/Pomeroy/Selig
% - Used for comparison to "known" results, AIAA-2012-2781
AF_name   = {'nlr7301'; '5412'};                                            % NACA airfoils to load [####(#)]
AF_load   = [1; 0];                                                         % Load [1] or create [0] airfoil
numPtsTot = {'100'; '100'};                                                 % Number of points per airfoil
AF_flip   = [1; 1];                                                         % Flip the airfoil vertically [1 = NO, -1 = YES]
AF_scale  = [1; 0.32];                                                      % Scale compared to first airfoil
AF_angle  = [0; 20];                                                        % Angle of rotation of airfoil [deg]
AF_offset = [0   0;                                                         % Offset from origin [x y]
             1-0.053 -0.026-0.0096];

% User-defined knowns
Vinf = 1;                                                                   % Freestream velocity []  (just leave this at 1)
AoA  = 0;                                                                   % Angle of attack [deg]

% Number of airfoils
numAF = length(AF_name);                                                    % Number of airfoils

% Plotting flags
flagPlot = [1;          % Airfoil with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison
            0;          % Airfoil streamlines
            0;          % Pressure coefficient contour
            0];         % Velocity difference contour plots (only for sanity check)

%% XFOIL - CREATE/LOAD AIRFOIL

% PPAR menu options
PPAR.P  = '4';                                                              % "Panel bunching parameter"
PPAR.T  = '1';                                                              % "TE/LE panel density ratios"
PPAR.R  = '1';                                                              % "Refined area/LE panel density ratio"
PPAR.XT = '1 1';                                                            % "Top side refined area x/c limits"
PPAR.XB = '1 1';                                                            % "Bottom side refined area x/c limits"

XB       = cell(numAF,1);                                                   % Initialize X boundary points
YB       = cell(numAF,1);                                                   % Initialize Y boundary points
numPtsAF = zeros(numAF,1);                                                  % Number of points
numPanAF = zeros(numAF,1);                                                  % Number of panels

for i = 1:1:numAF
    PPAR.N  = numPtsTot{i};                                                 % "Number of panel nodes"
    
    % Call XFOIL function to obtain the following:
    % - Airfoil coordinates
    % - Pressure coefficient along airfoil surface
    % - Lift, drag, and moment coefficients
    [xFoilResults,success] = XFOIL_N(AF_name{i},PPAR,AoA,AF_load(i));       % Call the XFOIL function
    if (success == 0)
        return;
    end
    
    % Separate out results from XFOIL function results   
    XB{i,1} = AF_scale(i)*xFoilResults.XB;                                  % Boundary point X-coordinate
    YB{i,1} = AF_flip(i)*AF_scale(i)*xFoilResults.YB;                       % Boundary point Y-coordinate
    
    % Rotate the airfoil
    v      = [XB{i} YB{i}];                                                 % Concatenate for rotation
    R      = [cosd(-AF_angle(i)) -sind(-AF_angle(i));                       % Rotation matrix based on each airfoil's rotation angle (AF_angle)
              sind(-AF_angle(i)) cosd(-AF_angle(i))];
    vo     = (R*v')';                                                       % Rotate the airfoil data
    x_rot  = vo(:,1);                                                       % Extract rotated X boundary points
    y_rot  = vo(:,2);                                                       % Extract rotated Y boundary points
    
    XB{i,1} = x_rot + AF_offset(i,1);                                       % Offset the X boundary points by user-defined placement (AF_offset)
    YB{i,1} = y_rot + AF_offset(i,2);                                       % Offset the Y boundary points by user-defined placement (AF_offset)
    
    numPtsAF(i,1) = length(XB{i});                                          % Number of points for each airfoil
    numPanAF(i,1) = length(XB{i})-1;                                        % Number of panels for each airfoil
end

%% PANEL CALCULATIONS

totPan   = sum(numPtsAF)-1;                                                 % Total number of panels (including false panels)
falsePan = cumsum(numPanAF) + (1:1:numAF)';                                 % Indices of false panels
falsePan = falsePan(1:end-1);                                               % Get rid of last entry in array

% Array of actual airfoil panel indices
panInd           = (1:1:totPan)';                                           % List of all panels
panInd(falsePan) = 0;                                                       % Set intermediate panels to a value of zero
panInd           = panInd(panInd ~= 0);                                     % Panel indices of the actual panels (not intermediate panels)
iInd             = panInd;                                                  % Rename for use in functions
jInd             = panInd;                                                  % Rename for use in functions

% Starting/ending indices for referencing each airfoil's indices in panInd
kS(1,1) = 1;                                                                % Set starting index for first airfoil
for k = 1:1:numAF                                                           % Loop over all airfoils
    if (k == 1)                                                             % If it's the first airfoil
        kE(k,1) = kS(k) + numPanAF(k) - 1;                                  % Ending index of the first airfoil
    else                                                                    % If it's not the first airfoil
        kS(k,1) = kS(k-1) + numPanAF(k-1);                                  % Starting index of the rest of the airfoils
        kE(k,1) = kS(k) + numPanAF(k) - 1;                                  % Ending index of the rest of the airfoils
    end
end

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points on each airfoil
% - Flips only airfoils that are incorrectly oriented
for i = 1:1:numAF                                                           % Loop over all airfoils
    edge = zeros(numPanAF(i),1);                                            % Initialize edge value array
    for j = 1:1:numPanAF(i)                                                 % Loop over all panels
        edge(j,1) = (XB{i}(j+1)-XB{i}(j))*(YB{i}(j+1)+YB{i}(j));            % Compute edge values
    end
    sumEdge = sum(edge);                                                    % Sum all edge values
    
    if (sumEdge < 0)                                                        % If panels are CCW
        XB{i} = flipud(XB{i});                                              % Flip the X-data array
        YB{i} = flipud(YB{i});                                              % Flip the Y-data array
    end
end

% Individual airfoils to save for later plotting
XB_AF = XB;                                                                 % X boundary points for each airfoil
YB_AF = YB;                                                                 % Y boundary points for each airfoil

% Consolidate all airfoils into one array
XBT = [];                                                                   % Initialize temporary array
YBT = [];                                                                   % Initialize temporary array
for i = 1:1:numAF                                                           % Loop over all airfoils
    XBT = [XBT; XB{i}];                                                     % Concatenate all airfoil X boundary points
    YBT = [YBT; YB{i}];                                                     % Concatenate all airfoil Y boundary points
end
XB = XBT;                                                                   % Overwrite X boundary points
YB = YBT;                                                                   % Overwrite Y boundary points

% Number of panels
numPtsTot = length(XB);                                                     % Number of boundary/control points
numPanTot = length(XB)-1;                                                   % Number of panels

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
XC    = zeros(numPanTot,1);                                                 % Initialize control point X-coordinate array
YC    = zeros(numPanTot,1);                                                 % Initialize control point Y-coordinate array
S     = zeros(numPanTot,1);                                                 % Intialize panel length array
phiD  = zeros(numPanTot,1);                                                 % Initialize panel orientation angle array [deg]

% Find geometric quantities of the airfoils
for i = 1:1:numPanTot                                                       % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = sqrt(dx^2 + dy^2);                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face) [deg]
    if (phiD(i) < 0)                                                        % If panel angles are negative [deg]
        phiD(i) = phiD(i) + 360;                                            % Make all panel angles positive [deg]
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make all panel angles between 0 and 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

% Exclude panels that are not on airfoils
beta_Real = beta(panInd);                                                   % Keep values that are part of actual airfoils
S_Real    = S(panInd);                                                      % Keep values that are part of actual airfoils

%% COMPUTE SOURCE AND VORTEX PANEL STRENGTHS - REFS [2,3,6,7,10]

% Number of panels
numPan = sum(numPanAF);                                                     % Number of real airfoil panels

% Geometric integral
[I,J] = COMPUTE_IJ_SPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd);               % Call COMPUTE_IJ_SPM function (Refs [2] and [3])
[K,L] = COMPUTE_KL_VPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd);               % Call COMPUTE_KL_VPM function (Refs [6] and [7])

% Populate A matrix
% - Simpler option: A = I + pi*eye(numPan,numPan);
A = zeros(numPan,numPan);                                                   % Initialize the A matrix
for i = 1:1:numPan                                                          % Loop over all i panels
    for j = 1:1:numPan                                                      % Loop over all j panels
        if (i == j)                                                         % If the panels are the same
            A(i,j) = pi;                                                    % Set A equal to pi
        else                                                                % If panels are not the same
            A(i,j) = I(i,j);                                                % Set A equal to I
        end
    end
end

% Right column of A matrix
for k = 1:1:numAF                                                           % Loop over all airfoils
    for i = 1:1:numPan                                                      % Loop over all i panels (rows)
        A(i,numPan+k) = -sum(K(i,kS(k):kE(k)));                             % Add gamma term to right-most column of A matrix
    end
end

% Populate b array
b = zeros(numPan,1);                                                        % Initialize the b array
for i = 1:1:numPan                                                          % Loop over all i panels (rows)
    b(i) = -Vinf*2*pi*cos(beta_Real(i));                                    % Compute RHS array
end

% ===== Enforce the Kutta Condition, REF [10] =====
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

% Bottom row of A matrix (Kutta condition)
for k = 1:1:numAF                                                           % Loop over all airfoils
    for j = 1:1:numPan                                                      % Loop over all j panels (columns)
        A(numPan+k,j) = J(kS(k),j) + J(kE(k),j);                            % Source panel terms for each Kutta condition equation in A matrix
    end
    A(numPan+k,numPan+k) = -sum(L(kS(k),:) + L(kE(k),:)) + 2*pi;            % Vortex panel terms for each Kutta condition equation in A matrix
end

% Last element of b array (Kutta condition)
for k = 1:1:numAF                                                           % Loop over all airfoils
    b(numPan+k) = -Vinf*2*pi*(sin(beta_Real(kS(k)))+sin(beta_Real(kE(k)))); % Set b array value from Kutta condition equation for each airfoil
end

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
% ===== Enforce the Kutta Condition =====

% Compute result array
resArr = A\b;                                                               % Solution array including lambda and gamma

% Separate lambda and gamma values from result array
lambda = resArr(1:end-numAF);                                               % Source strengths for each panel
gamma  = resArr(end-numAF+1:end);                                           % Vortex strength for each airfoil

gammaArr = zeros(numPan,1);                                                 % Initialize the gamma array 
for k = 1:1:numAF                                                           % Loop over every airfoil
    gammaArr(kS(k):kE(k)) = gamma(k);                                       % Set vortex strength for each panel on each airfoil
end

sumsumG = sum(gammaArr.*S_Real);                                            % Sum of all the vortex strengths (times panel lengths)

% Chord length: Between minimum X coordinate and maximum X coordinate
chord_v2 = max(XBT) - min(XBT);                                             % Difference between max X value and min X value

% Chord length: Sum of wing and flap chords
chord_v3 = 0;
for k = 1:1:numAF
    chord_v3 = chord_v3 + (max(XB_AF{k}) - min(XB_AF{k}));                  % Sum chord lengths of all airfoils
end

% Display values to Command Window
fprintf('Sum of G  : %g\n',sumsumG);                                        % Sum of vortex panel strengths (times panel lengths)
fprintf('KJ Lift   : %g\n',2*sumsumG);                                      % Kutta-Joukowski lift (chord length of 1)
fprintf('KJ Lift 2 : %g\n',2*sumsumG/chord_v2);                             % Kutta-Joukowski lift (chord length using chord_v2)
fprintf('KJ Lift 3 : %g\n',2*sumsumG/chord_v3);                             % Kutta-Joukowski lift (chord length using chord_v3)

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities on each panel
Vt = zeros(numPan,1);                                                       % Initialize tangential velocity
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient
for i = 1:1:numPan
    term1 = Vinf*sin(beta_Real(i));                                         % Freestream term
    term2 = sum(lambda.*J(i,:)')/(2*pi);                                    % Source panel terms
    term3 = gammaArr(i)/2;                                                  % i == j vortex panel term
    term4 = -(gammaArr(i)/(2*pi)).*sum(L(i,:));                             % i ~= j vortex panel terms
    
    Vt(i) = term1 + term2 + term3 + term4;                                  % Compute total tangential velocity on panel i
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient on panel i
end

%% COMPUTE STREAMLINES - REFS [4,8]

if (flagPlot(5) == 1 || flagPlot(6) == 1 || flagPlot(7) == 1)
    
    % Grid parameters
    nGridX = 150;                                                           % X-grid for streamlines and contours
    nGridY = 150;                                                           % Y-grid for streamlines and contours
    xVals  = [min(XB)-0.5 max(XB)+0.5];                                     % X-grid extents [min, max]
    yVals  = [min(YB)-0.3 max(YB)+0.3];                                     % Y-grid extents [min, max]

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
    Vx = zeros(nGridX,nGridY);                                              % Initialize X velocity matrix
    Vy = zeros(nGridX,nGridY);                                              % Initialize Y velocity matrix
    
    % Find which points are in/on the airfoil
    inon = zeros(nGridX,nGridY);                                            % Initialize in/on matrix with zeros
    for m = 1:1:nGridX                                                      % Loop over all X grid points
        for n = 1:1:nGridY                                                  % Loop over all Y grid points
            inMax = 0;                                                      % Set value to zero before checking all airfoils
            onMax = 0;                                                      % Set value to zero before checking all airfoils
            for k = 1:1:numAF                                               % Loop over all airfoils
                [in,on] = inpolygon(XX(m,n),YY(m,n),XB_AF{k},YB_AF{k});     % Call function to determine if point is in/on the airfoil
                inMax = max(inMax,in);                                      % Update value based on new information for current airfoil
                onMax = max(onMax,on);                                      % Update value based on new information for current airfoil
            end
            if (inMax == 1 || onMax == 1)                                   % If the point is in/on any of the airfoils
                inon(m,n) = 1;                                              % Set the value equal to one (otherwise it will stay zero)
            end
        end
    end
    
    % Solve for grid point X and Y velocities
    fprintf('Computing streamlines\n');                                     % Display status in Command Window
    for m = 1:1:nGridX                                                      % Loop over all X grid points
        fprintf('X grid index: %i/%i\n',m,nGridX);                          % Display current iteration in Command Window
        for n = 1:1:nGridY                                                  % Loop over all Y grid points
            XP = XX(m,n);                                                   % Current iteration's X grid point
            YP = YY(m,n);                                                   % Current iteration's Y grid point
            [Mx,My] = STREAMLINE_SPM_N(XP,YP,XB,YB,phi,S,numPan,jInd);      % Compute Mx and My geometric integrals, REF [4]
            [Nx,Ny] = STREAMLINE_VPM_N(XP,YP,XB,YB,phi,S,numPan,jInd);      % Compute Nx and Ny geometric integrals, REF [8]
            
            if (inon(m,n) == 0)                                             % If the grid point is not in/on any airfoil
                term1 = Vinf*cosd(AoA) + sum(lambda.*Mx/(2*pi));            % Uniform flow term and contribution from source panels
                term2 = 0;                                                  % Reset value to zero for every grid point
                for k = 1:1:numAF                                           % Loop over all airfoils
                    term2 = term2 + sum(-gamma(k).*Nx(kS(k):kE(k))/(2*pi)); % Contribution from vortex panels
                end
                Vx(m,n) = term1 + term2;                                    % Combine two terms to get X-direction velocity
                
                term3 = Vinf*sind(AoA) + sum(lambda.*My/(2*pi));            % Uniform flow term and contribution from source panels
                term4 = 0;                                                  % Reset value to zero for every grid point
                for k = 1:1:numAF                                           % Loop over all airfoils
                    term4 = term4 + sum(-gamma(k).*Ny(kS(k):kE(k))/(2*pi)); % Contribution from source panels
                end
                Vy(m,n) = term3 + term4;                                    % Combine two terms to get Y-direction velocity
            end
        end
    end
    
    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector []
    CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
end

%% CIRCULATION AND SOURCE/VORTEX STRENGTH CHECKS - REF [12]

if (flagPlot(5) == 1 || flagPlot(6) == 1 || flagPlot(7) == 1)
    
    xExtents = diff(xVals);                                                 % Difference between left and right edge (X) of streamline region
    yExtents = diff(yVals);                                                 % Difference between top and bottom edge (Y) of streamline region
    
    % Compute circulation
    aa   = xExtents/2.25;                                                   % Ellipse horizontal half-length
    bb   = yExtents/2.25;                                                   % Ellipse vertical half-length
    x0   = xVals(1) + xExtents/2;                                           % Ellipse center X-coordinate
    y0   = yVals(1) + yExtents/2;                                           % Ellipse center Y-coordinate
    numT = 5000;                                                            % Number of points on ellipse
    [Circulation,xC,yC,VxC,VyC] = COMPUTE_CIRCULATION(aa,bb,x0,y0,numT,...  % Compute circulation around ellipse, REF [11]
                                                        Vx,Vy,XX,YY);
    
    % Print values to Command Window
    fprintf('======= CIRCULATION RESULTS =======\n');
    fprintf('Sum of L     : %3.7f\n',sum(lambda.*S_Real));                  % Print sum of source strengths ([L]ambda)
    fprintf('Sum of G     : %3.7f\n',sumsumG);                              % Print sum of vortex strengths ([G]amma)
    fprintf('Circulation  : %3.7f\n',Circulation);                          % Print circulation from ellipse
    fprintf('K-J from G   : %3.7f\n',2*sumsumG/chord_v2);                   % K-J lift coefficient using G and chord_v2
    fprintf('K-J from Circ: %3.7f\n',2*Circulation/chord_v2);               % K-J lift coefficient using circulation and chord_v2
    fprintf('===================================\n');
end

%% PLOTTING

% FIGURE: Airfoil with panel normal vectors
if (flagPlot(1) == 1)
    S_Real = S(panInd);                                                     % Real panel lengths
    betaD2 = betaD(panInd);                                                 % Real panel angles
    XC2    = XC(panInd);                                                    % Real panel X control points
    YC2    = YC(panInd);                                                    % Real panel Y control points
    
    figure(1);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoil bodies as filled black polygons
    end
    for i = 1:1:numPan                                                      % Loop over all panels
        X(1) = XC2(i);                                                      % Set X start of panel orientation vector
        X(2) = XC2(i) + S_Real(i)*cosd(betaD2(i)+AoA);                      % Set X end of panel orientation vector
        Y(1) = YC2(i);                                                      % Set Y start of panel orientation vector
        Y(2) = YC2(i) + S_Real(i)*sind(betaD2(i)+AoA);                      % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',2);                                       % Plot panel normal vector
    end
    if (flagPlot(5) == 1 || flagPlot(6) == 1)                               % If we computed streamlines and circulation
        plot(xC,yC,'k-','LineWidth',2);                                     % Plot the ellipse that the circulation was computed around
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    axis equal;                                                             % Set axes equal
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Control and boundary points
if (flagPlot(2) == 1)
    figure(2);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XB,YB,'k-','LineWidth',2);                                         % Solid black line connecting boundary points
    pB = plot(XB,YB,'ko','MarkerFaceColor','k');                            % Boundary points as black-filled circles
    pC = plot(XC,YC,'ko','MarkerFaceColor','r');                            % Control points as red-filled circles
	xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Cp vectors at airfoil control points
if (flagPlot(3) == 1)
    betaD2 = betaD(panInd);                                                 % Real panel angles
    XC2    = XC(panInd);                                                    % Real panel X control points
    YC2    = YC(panInd);                                                    % Real panel Y control points
    sf     = 0.25;                                                          % Cp scale factor
    
    figure(3);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    Cps = abs(Cp*sf);                                                       % Scale and make positive all Cp values
    for i = 1:1:length(Cps)                                                 % Loop over all panels
        X(1) = XC2(i);                                                      % Control point X-coordinate
        X(2) = XC2(i) + Cps(i)*cosd(betaD2(i)+AoA);                         % Ending X-value based on Cp magnitude
        Y(1) = YC2(i);                                                      % Control point Y-coordinate
        Y(2) = YC2(i) + Cps(i)*sind(betaD2(i)+AoA);                         % Ending Y-value based on Cp magnitude
        
        if (Cp(i) < 0)                                                      % If pressure coefficient is negative
            p{1} = plot(X,Y,'r-','LineWidth',2);                            % Plot as a red line
        elseif (Cp(i) >= 0)                                                 % If pressure coefficient is zero or positive
            p{2} = plot(X,Y,'b-','LineWidth',2);                            % Plot as a blue line
        end
    end
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoil bodies
    end
    try                                                                     % Make sure there are both positive and negative Cp values
        legend([p{1},p{2}],{'Negative Cp','Positive Cp'});                  % Show legend
    catch                                                                   % If there only positive or negative Cp values, don't add legend
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient
if (flagPlot(4) == 1)
	XC2 = XC(panInd);                                                       % Real panel X control points
    
    figure(4);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XC2,Cp,'ks','MarkerFaceColor','k');                                % Plot Cp for all airfoils
    xlabel('X Coordinate');                                                 % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim('auto');                                                          % Set X-axis limits
    ylim('auto');                                                           % Set Y-axis limits to auto
    set(gca,'Ydir','reverse')                                               % Reverse direction of Y-axis
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Airfoil streamlines
if (flagPlot(5) == 1)
    figure(5);                                                              % Create the figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    for i = 1:1:length(Ysl)                                                 % Loop over all Y streamline starting points
        sl = streamline(XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % Plot the streamline
        set(sl,'LineWidth',2,'Color','b');                                  % Set streamline line width
    end
%     quiver(XX,YY,Vx,Vy,'r');                                                % Plot velocity vector arrows (comment back in)
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoil bodies
    end
	xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient contour
if (flagPlot(6) == 1)
    figure(6);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    contourf(XX,YY,CpXY,100,'EdgeColor','none');                            % Plot Cp contour
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoil bodies
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Symmetric airfoil velocity differences
% - Only used while I was checking that the symmetric airfoil cases were correct
% - Vx, Vy, and Cp should be the same in the top half and bottom half of the domain
% - Plotting differences between top and bottom half
% - Only works for 150 point domain
if (flagPlot(7) == 1)
    VxU       = Vx(1:75,:);                                                 % X-velocity in upper half of domain
    VxL       = flipud(Vx(76:end,:));                                       % X-velocity in lower half of domain
    VxD       = abs(abs(VxU) - abs(VxL));                                   % Difference of X-velocity half-domains
    VxDiff    = [VxD; flipud(VxD)];                                         % Concatenate difference half-domains
    VxDiff    = (VxDiff./(max(max(Vx))))*100;                               % Convert difference to percentage
    
    VyU       = Vy(1:75,:);                                                 % Y-velocity in upper half of domain
    VyL       = flipud(Vy(76:end,:));                                       % Y-velocity in lower half of domain
    VyD       = abs(abs(VyU) - abs(VyL));                                   % Difference of Y-velocity half-domains
    VyDiff    = [VyD; flipud(VyD)];                                         % Concatenate difference half-domains
    VyDiff    = (VyDiff./(max(max(Vy))))*100;                               % Convert difference to percentage
    
    CpU       = CpXY(1:75,:);                                               % Cp in upper half of domain
    CpL       = flipud(CpXY(76:end,:));                                     % Cp in lower half of domain
    CpD       = abs(abs(CpU) -abs(CpL));                                    % Difference of Cp half-domains
    CpDiff    = [CpD; flipud(CpD)];                                         % Concatenate difference half-domains
    CpDiff    = (CpDiff./(max(max(Cp))))*100;                               % Convert difference to percentage
    
    figure(7);                                                              % Create the figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    
    subplot(3,1,1);                                                         % Subplot 1
    cla; hold on;                                                           % Get ready for plotting
    contourf(XX,YY,VxDiff,100,'EdgeColor','none');                          % Plot Vx difference contour
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoils
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    title('Vx Difference');                                                 % Set title
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
    colorbar;                                                               % Turn on colorbar
    
    subplot(3,1,2);                                                         % Subplot 2
    cla; hold on;                                                           % Get ready for plotting
    contourf(XX,YY,VyDiff,100,'EdgeColor','none');                          % Plot Vy difference contour
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoils
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    title('Vy Difference');                                                 % Set title
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
    colorbar;                                                               % Turn on colorbar
    
    subplot(3,1,3);                                                         % Subplot 3
    cla; hold on;                                                           % Get ready for plotting
    contourf(XX,YY,CpDiff,100,'EdgeColor','none');                          % Plot Vy difference contour
    for i = 1:1:numAF                                                       % Loop over all airfoils
        fill(XB_AF{i},YB_AF{i},'k');                                        % Plot airfoils
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    title('Cp Difference');                                                 % Set title
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
    colorbar;                                                               % Turn on colorbar
end

%% PAPER COMPARISON

% arrAoA             = [0; 2.5; 5; 7.5; 10; 12.5];                            % Angle of attack [deg]
% arrPaper           = [1.68; 2; 2.3; 2.6; 2.875; 3.1];                       % Lift coefficient from paper (approximated)
% arr_5412_20_c1     = [2.16042; 2.51976; 2.87432; 3.2234; 3.56634; 3.9025];  % Lift coefficient based on chord of 1
% arr_5412_20_c1p247 = [1.73132; 2.0193; 2.30343; 2.58318; 2.85801; 3.1274];  % Lift coefficient based on chord of 1.247 (chord_v2)
% arr_5412_20_c1p301 = [1.66043; 1.9366; 2.20912; 2.47741; 2.74099; 2.9994];  % Lift coefficient based on chord of 1.301 (chord_v3)
% 
% figure(8);                                                                  % Create the figure
% cla; hold on; grid on;                                                      % Get ready for plotting
% set(gcf,'Color','White');                                                   % Set color to white
% set(gca,'FontSize',12);                                                     % Set font size
% p{1} = plot(arrAoA,arrPaper,'k-','LineWidth',2);                            % Plot paper data
% p{2} = plot(arrAoA,arr_5412_20_c1,'rs','MarkerFaceColor','r');              % Plot c = 1 data
% p{3} = plot(arrAoA,arr_5412_20_c1p247,'bs','MarkerFaceColor','b');          % Plot c = 1.247 data
% p{4} = plot(arrAoA,arr_5412_20_c1p301,'ms','MarkerFaceColor','m');          % Plot c = 1.301 data
% legend([p{:}],{'Paper','Chord = 1','Chord = 1.247','Chord = 1.301'},...     % Add legend
%                 'Location','NW');
% xlabel('Angle of Attack [deg]');                                            % X-axis label
% ylabel('Lift Coefficient');                                                 % Y-axis label
% xlim([-1 13]);                                                              % X-axis limits
% ylim('auto');                                                               % Y-axis limits
% box on;                                                                     % Add box around plot
% zoom reset;                                                                 % Reset the zoom







