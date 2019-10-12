% TESTING THE LOAD_AIRFOIL.M FUNCTION
% Written by: JoshTheEngineer
% Started   : 01/17/19
% Updated   : 01/17/19 - Started code
%                      - Works as expected

clear;
clc;

% Define filepath and all filenames in that path
flpth = '.\Airfoil_DAT\';
% flpth = '.\Airfoil_DAT_Selig\';
flnms = dir(flpth);

% Extract the filenames from the 'flnms' structure
for i = 1:1:length(flnms)-2
    flnmArr{i,1} = flnms(i+2).name(1:end-4);
end

% Loop through all files and find good/bad airfoils
goodFoil = [];
badFoil  = [];
for kk = 1:1:length(flnmArr)
    fprintf('Iteration: %i/%i\n',kk,length(flnmArr));
    fprintf('\tFile: %s\n',flnmArr{kk});
    
    try
        [dataX,dataY] = LOAD_AIRFOIL(flnmArr{kk});
%         [dataX,dataY] = LOAD_AIRFOIL_SELIG(flnmArr{kk});
        
        goodFoil = [goodFoil; flnmArr(kk)];
        
        % Plot the airfoil
        figure(1);
        cla; hold on; grid on;
        set(gcf,'Color','White');
        plot(dataX,dataY,'ko-');
        xlim([0 1]);
        ylim('auto');
        axis equal;
    catch
        badFoil = [badFoil; flnmArr(kk)];
    end
end