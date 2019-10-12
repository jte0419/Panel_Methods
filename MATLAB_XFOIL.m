% CALLING XFOIL FROM MATLAB
% Written by: JoshTheEngineer
% YouTube: www.youtube.com/joshtheengineer
% Website: www.joshtheengineer.com
% Started: 01/01/19
% Updated: 01/01/19 - Started code
%                   - Works as expected

clear;
clc;

NACA       = '2412';                                                        % NACA 4-digit airfoil [str]
AoA        = '0';                                                           % Angle of attack [deg]
numNodes   = '35';                                                          % Panel nodes [#]
saveFlnmAF = 'Save_Airfoil.txt';                                            % Airfoil coordinates filename
saveFlnmCp = 'Save_Cp.txt';                                                 % Pressure coefficient filename

% Delete files if they exist
if (exist(saveFlnmAF,'file'))                                               % If airfoil coordinates file exists
    delete(saveFlnmAF);                                                     % Delete this file
end
if (exist(saveFlnmCp,'file'))                                               % If Cp file exists
    delete(saveFlnmCp);                                                     % Delete this file
end

% Create the airfoil
fid = fopen('xfoil_input.txt','w');                                         % Open the XFoil input file
fprintf(fid,['NACA ' NACA '\n']);                                           % Specify NACA airfoil
fprintf(fid,'PPAR\n');                                                      % Enter paneling menu
fprintf(fid,['N ' numNodes '\n']);                                          % Specify number of panel nodes
fprintf(fid,'\n\n');                                                        % Set nodes and exit to XFOIL menu

% Save the airfoil data points
fprintf(fid,['PSAV ' saveFlnmAF '\n']);                                     % Save airfoil coordinates to test file

% Find the Cp vs. X plot
fprintf(fid,'OPER\n');                                                      % Enter OPER menu
fprintf(fid,['Alfa ' AoA '\n']);                                            % Set angle of attack [deg]
fprintf(fid,['CPWR ' saveFlnmCp]);                                          % Write Cp data to text file

% Close file
fclose(fid);                                                                % Close the file

% Run XFoil using input file
cmd = 'xfoil.exe < xfoil_input.txt';                                        % String to feed to system function
[status,result] = system(cmd);                                              % Run XFoil with input text file

%% READ DATA FILE: AIRFOIL

saveFlnmAF = 'Save_Airfoil.txt';                                            % File name
fidAirfoil = fopen(saveFlnmAF);                                             % Open file for reading

dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...              % Read data from file
                                 'Delimiter','','HeaderLines',0);
fclose(fidAirfoil);                                                         % Close the file
delete(saveFlnmAF);                                                         % Delete the file

% Separate boundary points
XB = dataBuffer{1}(:,1);                                                    % Boundary point X-coordinate
YB = dataBuffer{1}(:,2);                                                    % Boundary point Y-coordinate

%% READ DATA FILE: PRESSURE COEFFICIENT

saveFlnmCp = 'Save_Cp.txt';                                                 % File name
fidCP = fopen(saveFlnmCp);                                                  % Open file for reading
dataBuffer = textscan(fidCP,'%f %f %f','HeaderLines',3,...                  % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
fclose(fidCP);                                                              % Close the file
delete(saveFlnmCp);                                                         % Delete the file

% Separate Cp data
X_0  = dataBuffer{1,1}(:,1);                                                % X-coordinate
Y_0  = dataBuffer{1,1}(:,2);                                                % Y-coordinate
Cp_0 = dataBuffer{1,1}(:,3);                                                % Cp data

%% PLOT DATA

% Split airfoil into (U)pper and (L)ower
XB_U = XB(YB >= 0);                                                         % Upper airfoil X boundary points
XB_L = XB(YB < 0);                                                          % Lower airfoil X boundary points
YB_U = YB(YB >= 0);                                                         % Upper airfoil Y boundary points
YB_L = YB(YB < 0);                                                          % Lower airfoil Y boundary points

% Split Xfoil results into (U)pper and (L)ower
Cp_U = Cp_0(YB >= 0);                                                       % Upper airfoil Cp values
Cp_L = Cp_0(YB < 0);                                                        % Lower airfoil Cp values
X_U  = X_0(YB >= 0);                                                        % Upper airfoil X coordinate
X_L  = X_0(YB < 0);                                                         % Lower airfoil X coordinate

% Plot: Airfoil
figure(1);                                                                  % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
plot(XB_U,YB_U,'b.-');                                                      % Plot upper airfoil boundary points
plot(XB_L,YB_L,'r.-');                                                      % Plot lower airfoil boundary points
xlabel('X Coordinate');                                                     % Set X-label
ylabel('Y Coordinate');                                                     % Set Y-label
axis equal;                                                                 % Set axes equal

% Plot: Pressure coefficient
figure(2);                                                                  % Create figure
cla; hold on; grid on;                                                      % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
plot(X_U,Cp_U,'bo-','LineWidth',2);                                         % Plot upper surface Cp values
plot(X_L,Cp_L,'ro-','LineWidth',2);                                         % Plot lower surface Cp values
xlabel('X Coordinate');                                                     % Set X-label
ylabel('Cp');                                                               % Set Y-label
ylim('auto');                                                               % Set Y-limits
set(gca,'Ydir','reverse')                                                   % Reverse Y-axis
