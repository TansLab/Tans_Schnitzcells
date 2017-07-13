
%%
%load('D:\2017_data_analyses_Giulia\2017-05-17_asc824_mCherry_calibration\antibiotics_part\pos1crop\data\pos1crop-Schnitz');
load('D:\2017_data_analyses_Giulia\2017-06-06_asc825_mVenus_calibration\antibiotics\pos2crop\data\pos2crop-Schnitz');

%save('D:\2017_data_analyses_Giulia\2017-05-17_asc824_mCherry_calibration\antibiotics_part\pos1crop\data\pos1crop-Schnitz','schnitzcells');
%%
% figure(1); clf; hold on;
% plot(schnitzcells(100).time,schnitzcells(100).length_fitNew,'o');
% 
% figure(2); clf; hold on;
% plot(schnitzcells(100).R_time,schnitzcells(100).R6_mean,'o');
% 
% %%
% figure(3); clf; hold on;
% plot(schnitzcells(100).time_atdR,schnitzcells(100).dR5_sum,'o');

%%
allTimesFluor = [schnitzcells(:).Y_time]; % time points that correspond to taking fluor image
allFluor = [schnitzcells(:).Y6_mean];     % post-processed value of fluorescence mean value in cell

allFluorSum2 = [schnitzcells(:).Y2_sum];     % post-processed value of fluorescence mean value in cell

allTimesDelta = [schnitzcells(:).time_atdY]; % time point right in between two fluor images used for delta
allFluorDelta = [schnitzcells(:).dY5_sum];   % for 1 cell, fluor image of frame N+1 minus fluor image of frame N

allTimesGrowth = [schnitzcells(:).time]; % for 1 cell, all times at which phase image was taken
allGrowth = [schnitzcells(:).muP9_fitNew_all]; % for 1 cell, for each frame, growth rate based on exponential growth base 2

allLengths = [schnitzcells(:).length_fitNew]; %total legnth of cells at certain point in time
allE = [schnitzcells(:).length_fitNew];
%%
figure(1); clf; hold on;
f = fit(allTimesFluor(:),allFluor(:),'exp1')
plot(f,allTimesFluor(:),allFluor(:),'g*');
title('Mean fluorescence (a.u.) vs. time (min)');
xlabel('time (min)');
ylabel('Mean fluorescence (a.u.)');

%%
figure(2); clf; hold on;
plot(allTimesDelta,allFluorDelta,'b*');
title('Relative increse in fluorescence (a.u.) vs. time (min)');
xlabel('DeltaTime (min)');
ylabel('DeltaFluorescence (a.u.)');
grid on;

uniqueTimesDelta = unique(allTimesDelta);
myBins = [uniqueTimesDelta-(uniqueTimesDelta(2)-uniqueTimesDelta(1))/2 uniqueTimesDelta(end)-(uniqueTimesDelta(2)+uniqueTimesDelta(1))/2]
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesDelta},{real(allFluorDelta)},myBins)
plot(binCenters,meanValuesForBins,'r--','LineWidth',1.5)

%%
figure(3); clf; hold on;
plot(allTimesGrowth,allGrowth,'r*');
title('Growth rate (dbl/hr) vs. time (min)');
xlabel('time (min)');
ylabel('Growth rate (dbl/hr)');

uniqueTimesGrowth = unique(allTimesGrowth);
myBins = [uniqueTimesGrowth-(uniqueTimesGrowth(2)-uniqueTimesGrowth(1))/2 uniqueTimesGrowth(end)-(uniqueTimesGrowth(2)+uniqueTimesGrowth(1))/2]
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesGrowth},{real(allGrowth)},myBins)
plot(binCenters,meanValuesForBins,'o-','LineWidth',3)

ylim([0,3]);

%%
figure(4); clf; hold on;
plot(allTimesGrowth,allLengths,'o');
title('Cells length vs. time (min)');
xlabel('time (min)');
ylabel('Cell Length');


uniqueTimesGrowth = unique(allTimesGrowth);
myBins = [uniqueTimesGrowth-(uniqueTimesGrowth(2)-uniqueTimesGrowth(1))/2 uniqueTimesGrowth(end)-(uniqueTimesGrowth(2)+uniqueTimesGrowth(1))/2]
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesGrowth},{real(allLengths)},myBins)
plot(binCenters,meanValuesForBins,'o-','LineWidth',3)

%%
figure(5); clf; hold on;
plot(allTimesFluor,allFluorSum2,'g*');
title('Summed fluorescence (a.u.) vs. time (min)');
xlabel('time (min)');
ylabel('Summed fluorescence (a.u.)');
% 
% uniqueTimesFluor = unique(allTimesFluor);
% myBins = [uniqueTimesFluor-(uniqueTimesFluor(2)-uniqueTimesFluor(1))/2 uniqueTimesFluor(end)-(uniqueTimesFluor(2)+uniqueTimesFluor(1))/2]
% [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=...
%     binnedaveraging({allTimesFluor},{real(allFluorSum2)},myBins)
% plot(binCenters,meanValuesForBins,'o--','LineWidth',1.5)
%%

figure(6); clf; hold on;
SCHNITZNR=50;
plot(schnitzcells(SCHNITZNR).R_time,schnitzcells(SCHNITZNR).R2_sum,'o-');
SCHNITZNR=100;
plot(schnitzcells(SCHNITZNR).R_time,schnitzcells(SCHNITZNR).R2_sum,'o-');
SCHNITZNR=30;
plot(schnitzcells(SCHNITZNR).R_time,schnitzcells(SCHNITZNR).R2_sum,'o-');

% 
timeOfSwitch   =datenum(2017,05,19,15,00,40);
time0InSchnitz =min([schnitzcells(:).timestamp]);

timeOfSwitchRelative = timeOfSwitch-time0InSchnitz;
timeOfSwitchRelativeMins = timeOfSwitchRelative*24*60;

figure(2);
plot([timeOfSwitchRelativeMins timeOfSwitchRelativeMins],[min(allFluorDelta),max(allFluorDelta)],':k')
%%
figure(3);
plot([timeOfSwitchRelativeMins timeOfSwitchRelativeMins],[0,3],':k','LineWidth',2)



