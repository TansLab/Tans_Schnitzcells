
%% 1A
% Some code to analyze growth curve manually
% %%%
% Note that some data needs to be pasted from Excel sheet w. Matlab
% commands!

% Set position here!
mypos = 'pos9crop';

close all;

% Load the following file to get the Schnitzes:
load(['D:\MICROSCOPE_EXPERIMENTS\To_Analyze\2014-04-04\' mypos '\data\' mypos '-Schnitz.mat'])
    % Todo: check difference with pos7crop_lin.mat (created 1 min earlier)
    
% prepare data
p1 = DJK_initschnitz(mypos,'2014-04-04','e.coli.AMOLF','rootDir','D:\MICROSCOPE_EXPERIMENTS\To_Analyze\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none','badseglistname','badseglist');
[p1,schnitzcells] = DJK_compileSchnitzImproved_3colors(p1,'quickMode',1);

% make fit
fitTime = DJK_analyzeMu(p1, schnitzcells, 'xlim', [0 1000], 'onScreen', 0,'fitTime',[0 200], 'dumpPlot', 1);
    % Last parameter, dumpPlot, is a custom MW one, which creates
    % parameters in the workspace with the plot.
    
% Plot the data
semilogy(data_frames, data_muField_sum,'o')









%%










