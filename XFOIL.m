function [xFoilResults,success] = XFOIL(NACA,PPAR,AoA,flagAirfoil)

% PURPOSE
% - Create or load airfoil based on flagAirfoil
% - Save and read airfoil coordinates
% - Save and read airfoil pressure coefficient
% - Save and read airfoil lift, drag, and moment coefficients
%
% INPUTS
% - NACA         : Four-digit NACA airfoil designation
% - PPAR         : Paneling variables used in XFOIL PPAR menu
% - AoA          : Angle of attack [deg]
% - flagAirfoil  : Flag for loading/creating airfoil
% 
% OUTPUTS
% - xFoilResults : Structure containing all results
% - success      : Flag indicating whether solution was a success

%% CALL XFOIL FROM MATLAB

xFoilResults = [];                                                          % Initialize results structure

if (flagAirfoil.XFoilCreate == 1)                                           % If the user wants XFOIL to create a NACA airfoil
    airfoilName         = NACA;                                             % Set the airfoilName to the input NACA airfoil
    xFoilResults.afName = airfoilName;                                      % Send the airfoil name back from this function
    success             = 1;                                                % This will be successful
elseif (flagAirfoil.XFoilLoad == 1)                                         % If the user wnats to load a DAT file airfoil
    [flnm,~,success]    = uigetfile('./Airfoil_DAT_Selig/*.dat',...         % User input of airfoil file to load
                                    'Select Airfoil File');
    airfoilName         = flnm(1:end-4);                                    % Set the airfoilName based on loaded file
    xFoilResults.afName = airfoilName;                                      % Send the airfoil name back from this function
    if (success == 0)                                                       % If the user exited dialog box without selecting airfoil
        return;                                                             % Exit the function
    else                                                                    % If user selected an airfoil
        success = 1;                                                        % This will be successful
    end
end

% Save-to file names
saveFlnm    = ['Save_' airfoilName '.txt'];                                 % Airfoil coordinates save-to file
saveFlnmCp  = ['Save_' airfoilName '_Cp.txt'];                              % Airfoil Cp save-to file
saveFlnmPol = ['Save_' airfoilName '_Pol.txt'];                             % Airfoil polar save-to file

% Delete files if they exist
if (exist(saveFlnm,'file'))                                                 % If airfoil coordinate file exists
    delete(saveFlnm);                                                       % Delete it
end
if (exist(saveFlnmCp,'file'))                                               % If airfoil Cp file exists
    delete(saveFlnmCp);                                                     % Delete it
end
if (exist(saveFlnmPol,'file'))                                              % If airfoil polar file exists
    delete(saveFlnmPol);                                                    % Delete it
end

% Create the airfoil
fid = fopen('xfoil_input.inp','w');                                         % Create an XFoil input file, and make it write-able
if (flagAirfoil.XFoilLoad == 1)                                             % If user wants to load DAT airfoil file
    fprintf(fid,['LOAD ' './Airfoil_DAT_Selig/' flnm '\n']);                % Load selected airfoil
elseif (flagAirfoil.XFoilCreate == 1)                                       % If user wants to specify a 4-digit airfoil
    fprintf(fid,['NACA ' NACA '\n']);                                       % Specify NACA airfoil
end
fprintf(fid,'PPAR\n');                                                      % Enter the PPAR (paneling) menu
fprintf(fid,['N ' PPAR.N '\n']);                                            % Define "Number of panel nodes"
fprintf(fid,['P ' PPAR.P '\n']);                                            % Define "Panel bunching paramter"
fprintf(fid,['T ' PPAR.T '\n']);                                            % Define "TE/LE panel density ratios"
fprintf(fid,['R ' PPAR.R '\n']);                                            % Define "Refined area/LE panel density ratio"
fprintf(fid,['XT ' PPAR.XT '\n']);                                          % Define "Top side refined area x/c limits"
fprintf(fid,['XB ' PPAR.XB '\n']);                                          % Define "Bottom side refined area x/c limits"
fprintf(fid,'\n');                                                          % Apply all changes
fprintf(fid,'\n');                                                          % Back out to XFOIL menu

% Save the airfoil data points
fprintf(fid,['PSAV ' saveFlnm '\n']);                                       % Save the airfoil coordinate file

% Get Cp and polar data
fprintf(fid,'OPER\n');                                                      % Enter OPER menu
fprintf(fid,'Pacc 1 \n');                                                   % Begin polar accumulation
fprintf(fid,'\n\n');                                                        % Don't enter save or dump file names
fprintf(fid,['Alfa ' num2str(AoA) '\n']);                                   % Set angle of attack
fprintf(fid,['CPWR ' saveFlnmCp '\n']);                                     % Write the Cp file
fprintf(fid,'PWRT\n');                                                      % Save the polar data
fprintf(fid,[saveFlnmPol '\n']);                                            % Save polar data to this file
if (exist(saveFlnmPol) ~= 0)                                                % If saveFlnmPol already exists
    fprintf(fid,'y \n');                                                    % Overwrite existing file
end

fclose(fid);                                                                % Close the input file

cmd = 'xfoil.exe < xfoil_input.inp';                                        % Write command to run XFoil

[~,~] = system(cmd);                                                        % Run XFoil with the input file

fclose all;                                                                 % Close all files

% Delete the XFoil input file
if (exist('xfoil_input.inp','file'))                                        % If the input file exists
    delete('xfoil_input.inp');                                              % Delete the file
end

%% READ CP DATA

fidCP = fopen(saveFlnmCp);                                                  % Open the airfoil Cp file
dataBuffer = textscan(fidCP,'%f %f %f','HeaderLines',3,...                  % Read in X, Y, and Cp data
                                     'CollectOutput',1,...
                                     'Delimiter','');
fclose(fidCP);                                                              % Close the file
delete(saveFlnmCp);                                                         % Delete the file

% Save airfoil Cp data to function solution variable
xFoilResults.X  = dataBuffer{1,1}(:,1);                                     % X-data points
xFoilResults.Y  = dataBuffer{1,1}(:,2);                                     % Y-data points
xFoilResults.CP = dataBuffer{1,1}(:,3);                                     % Cp data

%% READ AIRFOIL COORDINATES

fidAirfoil = fopen(saveFlnm);                                               % Open the airfoil file

dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...              % Read the XB and YB data from the data file
                                  'Delimiter','','HeaderLines',0);

XB = dataBuffer{1}(:,1);                                                    % Boundary point X-coordinate
YB = dataBuffer{1}(:,2);                                                    % Boundary point Y-coordinate
fclose(fidAirfoil);                                                         % Close the airfoil file
delete(saveFlnm);                                                           % Delete the airfoil file

% Save airfoil boundary points to function solution variable
xFoilResults.XB = XB;                                                       % Airfoil boundary X-points
xFoilResults.YB = YB;                                                       % Airfoil boundary Y-points

%% READ POLAR DATA

% Load and read polar file (Save_Polar.txt)
fidPolar = fopen(saveFlnmPol);                                              % Open polar data file for reading

dataBuffer = textscan(fidPolar,'%f %f %f %f %f %f %f','CollectOutput',1,... % Read data from file
                                 'Delimiter','','HeaderLines',12);
fclose(fidPolar);                                                           % Close the file
delete(saveFlnmPol);                                                        % Delete the file

% Extract polar data from buffer into function solution variable
xFoilResults.CL = dataBuffer{1,1}(2);                                       % Lift coefficient
xFoilResults.CD = dataBuffer{1,1}(3);                                       % Drag coefficient
xFoilResults.CM = dataBuffer{1,1}(5);                                       % Moment coefficient

