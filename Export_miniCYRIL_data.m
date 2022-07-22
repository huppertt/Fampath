% export miniCYRIL data

clear all
%% Load data

[file, path] = uigetfile('D:/Data/Data Analysis/miniCYRIL/miniCYRIL at Columbia/Processed NIRS Data/*.mat','Select processed file');
load([path file]);

%% create MATLAB table
VarNames = {'Timestamp','Time_s','HbO2_mM','HHb_mM','oxCCO_mM','Events'};
T = table(Timestamp, Time, Conc(:,1),Conc(:,1),Conc(:,1), Events(:,1), 'VariableNames',VarNames);

%% save table
filename = ['D:/Data/Data Analysis/miniCYRIL/miniCYRIL at Columbia/Exported NIRS Data/' erase(file,'.mat') ' table.mat'];
save(filename,'T')

%% export to .csv file

filename = ['D:/Data/Data Analysis/miniCYRIL/miniCYRIL at Columbia/Exported NIRS Data/' erase(file,' processed.mat') '.csv'];
writetable(T, filename);

