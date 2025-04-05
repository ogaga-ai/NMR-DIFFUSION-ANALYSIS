clear; clc; clf; close all;
set(0,'DefaultFigureVisible','on')
format shortg
warning off

%% Experiment Details
for z = 1
    [file, path] = uigetfile('*', 'Select Data File');
    filename = fullfile(path, file);
    data =  xlsread(filename);
    file
    
    chart_title = 'Input';  %Title outside dialog box
    dims = [1 35]; %Dimension of dialog box
    
    prompt = "# of Peaks";
    numPeakss = inputdlg(prompt,chart_title,dims);
    numPeaks = str2double(numPeakss)
    if isempty(numPeaks)
        fprintf('Code terminated: Missing Input.')
        return
    end
    
    prompt = "# of B Values (slices)";
    numPointss = inputdlg(prompt,chart_title,dims);
    numPoints = str2double(numPointss)
    if isempty(numPoints)
        fprintf('Code terminated: Missing Input.')
        return
    end
    
    P = data(:,3);
    skip = 22;
end

%% Extracting relevant Indices
for z = 1
    t = 100;
    n = mod(0:t-1,2);
    s = mod(1:t,2);
    idx(1) = 1;
    
    for i = 1:numPeaks*2
        idx(i+1) = idx(i) + s(i)*(numPoints-1) + n(i)*(skip);
    end
end

%% Organizing and Normalizing Data
for z = 1
    S = [];
    j = 1;
    for i = 1:2:numPeaks*2
        S(:,j) = P(idx(i):idx(i+1))/P(idx(i)); j = j+1;
    end
end

%% Print
for z = 1
    [~, name, ~] = fileparts(file);
    T_file_csv = [sprintf('Normalized_%s', name), '.csv'];
    T_name = fullfile(path, T_file_csv);
    writematrix(S, T_name)
    winopen(T_name)
end
