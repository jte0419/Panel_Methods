% FUNCTION - LOAD AIRFOIL
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 01/19/17
% Updated   : 01/19/17 - Started code
%                      - Works as expected for all but naca23015.dat (but
%                        that is expected)
%                      - Panel method code works with all files loaded
%                        using this function
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

function [dataX,dataY] = LOAD_AIRFOIL(fileName)

% Get the number of header lines in airfoil file
fidAirfoil = fopen(['./Airfoil_DAT/' fileName '.dat']);                     % Open the airfoil data file
dataCell   = textscan(fidAirfoil,'%s','Delimiter','\n');                    % Read each line of the file into a cell array of strings
fclose(fidAirfoil);                                                         % Close the airfoil data file
hdrlns     = nan;                                                           % Set the hdrlns variable to nan
for k = 1:1:length(dataCell{1})                                             % Loop over all lines of the file
    charArr = cell2mat(dataCell{1,1}(k));                                   % Conver the cell into a character array
    try
        C = strsplit(charArr,' ');                                          % Try to split the character array at any spaces
        if (str2double(C(1)) <= 1 && str2double(C(1)) >= 0)                 % If the first value of the split is between 0 and 1
            fprintf('\tHeader Lines: %i\n',k-1);                            % Then we have found the first data point of the airfoil
            hdrlns = k-1;                                                   % Set the number of header lines
            break;                                                          % Break out of the for loop
        end
    catch
    end
    if (isnan(hdrlns))                                                      % If we couldn't find any header lines
        C = strsplit(charArr);                                              % Split the character array
        if (str2double(C(1)) <= 1 && str2double(C(1)) >= 0)                 % Check if first character is a '1' or a '0'
            fprintf('\tHeader Lines: %i\n',k-1);                            % Then we have found the first data point of the airfoil
            hdrlns = k-1;                                                   % Set the number of header lines
            break;                                                          % Break out of the for loop
        end
    end
end

% Check whether we will have two data columns
firstLine = dataCell{1}(hdrlns+1);                                          % Get first line of actual data
split     = strsplit(cell2mat(firstLine),' ');                              % Split data at a space (should be between X and Y data columns)

% Load the data in, getting rid of empty lines
if (length(split) == 1)                                                     % If we don't have normal X and Y split data (1 column)
    dataX = zeros(length(dataCell{1})-hdrlns,1);                            % Initialize X data array
    dataY = zeros(length(dataCell{1})-hdrlns,1);                            % Initialize Y data array
    
    % Loop through the dataCell instead of parsing with textscan
    for i = 1:1:length(dataCell{1})-hdrlns                                  % Loop through dataCell
        if (length(dataCell{1}{i+hdrlns}) == 3)                             % If there are only three characters
            vals       = cell2mat(dataCell{1}(i+hdrlns));                   % Convert to double
            dataX(i,1) = str2double(vals(1));                               % Set the airfoil X-data
            dataY(i,1) = str2double(vals(3));                               % Set the airfoil Y-data
        elseif (length(dataCell{1}{i+hdrlns}) == 0)                         %#ok<*ISMT>
            % Do nothing
        else
            vals       = strsplit(dataCell{1}{i+hdrlns});                   % Split the string
            dataX(i,1) = str2double(vals(1));                               % Set the airfoil X-data
            dataY(i,1) = str2double(vals(2));                               % Set the airfoil Y-data
        end
    end
else                                                                        % If we have normal X and Y split data (2 columns)
    fidAirfoil = fopen(['./Airfoil_DAT/' fileName '.dat']);                 % Open the airfoil data file
    dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...          % Read in the airfoil data
                                             'Delimiter','',...
                                             'Whitespace',' ',...
                                             'HeaderLines',hdrlns);
    dataX = dataBuffer{1}(:,1);                                             % Airfoil X-data
    dataY = dataBuffer{1}(:,2);                                             % Airfoil Y-data
    fclose(fidAirfoil);                                                     % Close the airfoil data file
end

% Specific exception for 'bacj' airfoil
if (strcmpi('bacj',fileName))
    dataX = dataX(1:end-1);                                                 % Get rid of last data point
    dataY = dataY(1:end-1);                                                 % Get rid of last data point
end

% Find the two maximum values in the array
indMax = find(dataX == max(dataX));                                         % Will get the index (or indices) of maximum dataX values
if (length(indMax) < 2)                                                     % If we don't have 2 maximum value indices yet
    dataX2            = dataX;                                              % Set intermediate variable
    dataX2(indMax(1)) = 0;                                                  % Set the maximum value to '0'
    secInd            = find(dataX2 == max(dataX2));                        % Find the second maximum value index
    if (min(indMax > secInd) == 1)                                          % If indMax values are all bigger than secInd values
        indMax(2) = secInd(1);                                              % Set second maximum index to the first value in second index
    else                                                                    % Otherwise
        indMax(2) = secInd(secInd > indMax(1));                             % Set second maximum index to the index in secInd that's bigger than indMax
    end
end

% Specific exception for three airfoils
if (strcmp(fileName,'ste87151') || strcmp(fileName,'ste87391') || strcmp(fileName,'stf86361'))
    indMax = [indMax(1); length(dataX)];                                    % Set maximum indices to first and last index
end

% Make sure first index is smaller than second index
indMax = sort(indMax);

% Orient the data correctly
if (indMax(1) == 1 && indMax(2) == length(dataX))                           % Do nothing, already oriented correctly
    dataXX = dataX;
    dataYY = dataY;
elseif (indMax(1) == 1 && indMax(2) ~= length(dataX))                       % Flip second half
    dataXX = [dataX(1:indMax(2)-1); flipud(dataX(indMax(2):end))];
    dataYY = [dataY(1:indMax(2)-1); flipud(dataY(indMax(2):end))];
elseif (indMax(1) ~= 1 && indMax(2) == length(dataX))                       % Flip first half
    dataXX = [flipud(dataX(1:indMax(1))); dataX(indMax(1)+1:end)];
    dataYY = [flipud(dataY(1:indMax(1))); dataY(indMax(1)+1:end)];
elseif (indMax(1) ~= 1 && indMax(2) ~= length(dataX))                       % Flip both
    dataXX = [dataX(indMax(2):end); dataX(1:indMax(1))];
    dataYY = [dataY(indMax(2):end); dataX(1:indMax(1))];
end
dataX   = dataXX;                                                           % Reset X-data to correctly oriented data
dataY   = dataYY;                                                           % Reset Y-data to correctly oriented data
dataArr = [dataX dataY];                                                    % Create [Nx2] array of X and Y data

% Delete any duplicate (0,0) lines (only need one)
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
