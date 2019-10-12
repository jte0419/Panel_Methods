% FUNCTION - LOAD AIRFOIL (SELIG FORMAT)
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 01/23/19
% Updated   : 01/23/19 - Started code
%                      - Works as expected
% 
% Purpose
% - Load the X and Y coordinates for airfoil from fileName
% - Data comes from files on UIUC website
% - Loads every airfoil EXCEPT naca23015.dat (has weird data)
% - Checked to make sure the loaded data works with panel method codes
% 
% INPUTS
% - fileName : File name of airfoil without extension (must be .dat)
% 
% OUTPUTS
% - dataX : Airfoil X-coordinate data array [Nx1]
% - dataY : Airfoil Y-coordinate data array [Nx1]

function [dataX,dataY] = LOAD_AIRFOIL_SELIG(fileName)

hdrlns = 1;
if (strcmp(fileName,'nasasc2-0714'))
    hdrlns = 3;
elseif (strcmp(fileName,'s1020'))
    hdrlns = 2;
end

fidAirfoil = fopen(['./Airfoil_DAT_Selig/' fileName '.dat']);               % Open the airfoil data file
dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...              % Read in the airfoil data
                                         'HeaderLines',hdrlns,...
                                         'Delimiter','');
dataX = dataBuffer{1}(:,1);                                                 % Airfoil X-data
dataY = dataBuffer{1}(:,2);                                                 % Airfoil Y-data
fclose(fidAirfoil);                                                         % Close the airfoil data file

% Delete any duplicate (0,0) lines (only need one)
dataArr     = [dataX dataY];
[~,ia,~]    = unique(dataArr,'rows','stable');                              % Find the unique values of the array
i           = true(size(dataArr,1),1);                                      % Set every index to true
i(ia)       = false;                                                        % Set indices false that are not duplicates
dataArr(i,:)= [];                                                           % Get rid of duplicate rows (true in i)
dataX       = dataArr(:,1);                                                 % Separate out X-data
dataY       = dataArr(:,2);                                                 % Separate out Y-data

% Close the airfoil by adding final point
if (dataY(1) ~= dataY(end))                                                 % If the start and end points are not the same
    dataX(end+1) = dataX(1);                                                % Create a new X-endpoint to close the airfoil
    dataY(end+1) = dataY(1);                                                % Create a new Y-endpoint to close the airfoil
end
