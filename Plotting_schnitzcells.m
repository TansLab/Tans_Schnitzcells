
%%
%load('D:\2017_data_analyses_Giulia\2017-05-17_asc824_mCherry_calibration\antibiotics_part\pos1crop\data\pos1crop-Schnitz');
% load('D:\2017_data_analyses_Giulia\2017-05-17_asc824_mCherry_calibration\antibiotics_part\pos1crop\data\pos1crop-Schnitz');
load('D:\2017_data_analyses_Giulia\2017-06-28_asc976_rRNA experiment\TY experiment\pos1crop\data\pos1crop-Schnitz');

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
allTimesFluor = [schnitzcells(:).C_time]; % time points that correspond to taking fluor image
allFluor = [schnitzcells(:).C6_mean];     % post-processed value of fluorescence mean value in cell

allFluorSum2 = [schnitzcells(:).C2_sum];     % post-processed value of fluorescence summed all pixels value in cell
allFluorMean2 = [schnitzcells(:).C2_mean];     % post-processed value of fluorescence mean value in cell
% allFluorDelta2 = [schnitzcells(:).dY2]; 

allTimesDelta = [schnitzcells(:).time_atdC]; % time point right in between two fluor images used for delta
allFluorDelta = [schnitzcells(:).dC5_sum];   % for 1 cell, fluor image of frame N+1 minus fluor image of frame N


allTimesGrowth = [schnitzcells(:).time]; % for 1 cell, all times at which phase image was taken
allGrowth = [schnitzcells(:).muP9_fitNew_all]; % for 1 cell, for each frame, growth rate based on exponential growth base 2
 

allLengths = [schnitzcells(:).length_fitNew]; %total legnth of cells at certain point in time
allE = [schnitzcells(:).length_fitNew];

%variables for scatter plot
allFluorRel = [schnitzcells(:).dC5_cycCor];
allConc = [schnitzcells(:).C6_mean_cycCor];
allGrowthDeltaFluo = [schnitzcells(:).muP9_fitNew_atdC5_cycCor];
allGrowthFluo = [schnitzcells(:).muP9_fitNew_cycCor];
%%
figure(1); clf; hold on;
% plot(f,allTimesFluorx,allFluorx,'g*');
plot(allTimesFluor,allFluor,'b*','Markersize',3);
title('Mean fluorescence vs. time (ASC976 mCerulean)');
xlabel('time (min)');
ylabel('Mean fluorescence (a.u.)');
grid on
MW_makeplotlookbetter(20);

% t = ~isnan(allTimesFluor) & ~isnan(allFluor); %use in case of NaN values
% allTimesFluorx = allTimesFluor(:);
% allFluorx = allFluor(:);
% f = fit(allTimesFluorx(t),allFluorx(t),'exp1');
% f = fit(allTimesFluor(:),allFluor(:),'exp1')

uniqueTimesFluor = unique(allTimesFluor);

myBins = [uniqueTimesFluor-(uniqueTimesFluor(2)-uniqueTimesFluor(1))/2, uniqueTimesFluor(end)-(uniqueTimesFluor(2)-uniqueTimesFluor(1))/2]
[meanValuesForBins,medianValuesForBins,binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesFluor},{real(allFluor)},myBins)
plot(binCenters,meanValuesForBins,'r-.','LineWidth',1)
% plot(binCenters,medianValuesForBins,'b-','LineWidth',2, 'Markersize', 1)

% ylim([200 900]);
xlim([0 80]);
% xlim([0 80])
% ylim([100 900]);
% % % xlim([0 250]);
% xlim([0 360])
% line([120 120], [-80000000 80000000])
% line([240 240], [-800000000 800000000])


%%
figure(2); clf; hold on;
plot(allTimesDelta,allFluorDelta,'b*');
title('Deltafluo vs. time (ASC976 mCerulean)');
xlabel('DeltaTime (min)');
ylabel('DeltaFluorescence (a.u.)');
grid on;
MW_makeplotlookbetter(20);

uniqueTimesDelta = unique(allTimesDelta);
myBins = [uniqueTimesDelta-(uniqueTimesDelta(2)-uniqueTimesDelta(1))/2 uniqueTimesDelta(end)-(uniqueTimesDelta(2)+uniqueTimesDelta(1))/2]
[meanValuesForBins,medianValuesForBins,binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesDelta},{real(allFluorDelta)},myBins)
plot(binCenters,meanValuesForBins,'r-.','LineWidth',1)
% plot(binCenters,medianValuesForBins,'g-','LineWidth',2, 'Markersize', 1)

% line([25 25], [-800000 800000]) % draw line for antibiotics experiments
ylim([-10000 20000]);
xlim([0 80]);
% xlim([0 360])
% line([120 120], [-80000000 80000000])
% line([240 240], [-8000000 8000000])
%%
figure(3); clf; hold on;
plot(allTimesGrowth,allGrowth,'r*');
title('Growth rate (ASC976 mCerulean)');
xlabel('time (min)');
ylabel('Growth rate (dbl/hr)');
grid on
MW_makeplotlookbetter(20);
uniqueTimesGrowth = unique(allTimesGrowth);

myBins = [uniqueTimesGrowth-(uniqueTimesGrowth(2)-uniqueTimesGrowth(1))/2 uniqueTimesGrowth(end)-(uniqueTimesGrowth(2)+uniqueTimesGrowth(1))/2]
[meanValuesForBins,medianValuesForBins,binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesGrowth},{real(allGrowth)},myBins)

plot(binCenters,meanValuesForBins,'k-','LineWidth',2, 'Markersize', 1)
plot(binCenters,medianValuesForBins,'b-','LineWidth',2, 'Markersize', 1)
% 
ylim([-1 3.5]);
% % xlim([0 250]);
xlim([0 80])
% % line([120 120], [-80 80])
% line([240 240], [-80 80])
%%
figure(4); clf; hold on;
plot(allTimesGrowth,allLengths,'o');
title('Cells length vs. time (min)');
xlabel('time (min)');
ylabel('Cell Length');
grid on
MW_makeplotlookbetter(20);

uniqueTimesGrowth = unique(allTimesGrowth);
myBins = [uniqueTimesGrowth-(uniqueTimesGrowth(2)-uniqueTimesGrowth(1))/2 uniqueTimesGrowth(end)-(uniqueTimesGrowth(2)+uniqueTimesGrowth(1))/2]
[meanValuesForBins,medianValuesForBins,binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesGrowth},{real(allLengths)},myBins)

plot(binCenters,meanValuesForBins,'o--','LineWidth',1, 'Markersize', 1)
% plot(binCenters,medianValuesForBins,'R-','LineWidth',2, 'Markersize', 1)

ylim([1 4.5]);
% xlim([0 250]);
xlim([0 80])
% line([25 25], [-80 80])

% xlim([0 250]);
%%
figure(5); clf; hold on;
plot(allTimesFluor,allFluorSum2,'g*');
title('Summed fluorescence (a.u.) vs. time (min)');
xlabel('time (min)');
ylabel('Summed fluorescence (a.u.)');
grid on
MW_makeplotlookbetter(20);
uniqueTimesFluor = unique(allTimesFluor);
myBins = [uniqueTimesFluor-(uniqueTimesFluor(2)-uniqueTimesFluor(1))/2 uniqueTimesFluor(end)-(uniqueTimesFluor(2)+uniqueTimesFluor(1))/2]
[meanValuesForBins,medianValuesForBins,binCenters,stdValuesForBins,stdErrValuesForBins]=...
    binnedaveraging({allTimesFluor},{real(allFluorSum2)},myBins)

plot(binCenters,meanValuesForBins,'o--','LineWidth',1, 'Markersize', 1)
% plot(binCenters,medianValuesForBins,'b-','LineWidth',2, 'Markersize', 1)

% ylim([-10000 20000]);
% xlim([0 250]);
% xlim([0 80])
% line([25 25], [-200000 1000000])
%%

figure(6); clf; hold on;


%SCHNITZNR=50; SCHNITZNR=100; SCHNITZNR=30;
SCHNITZNRS=[50, 100, 30];

SCHNITZNRS=1:numel(schnitzcells);

someColors=linspecer(numel(SCHNITZNRS));

for idx=1:numel(SCHNITZNRS)
    currentSchnitz=SCHNITZNRS(idx)
    plot(schnitzcells(currentSchnitz).time_atC,schnitzcells(currentSchnitz).C5_sum,...
        'o-','Color',someColors(idx,:),'LineWidth',2);
    plot(schnitzcells(currentSchnitz).time_atdC,schnitzcells(currentSchnitz).dC5_sum,...
        's-','Color',someColors(idx,:),'LineWidth',2);
end

plot([0,80],[0,0],'k-','LineWidth',2)
%%
timeOfSwitch   =datenum(2017,05,19,15,00,40);
time0InSchnitz =min([schnitzcells(:).timestamp]);

timeOfSwitchRelative = timeOfSwitch-time0InSchnitz;
timeOfSwitchRelativeMins = timeOfSwitchRelative*24*60;

figure(2);
plot([timeOfSwitchRelativeMins timeOfSwitchRelativeMins],[min(allFluorDelta),max(allFluorDelta)],':k')
%%
figure(3);
plot([timeOfSwitchRelativeMins timeOfSwitchRelativeMins],[0,3],':k','LineWidth',2)

%%
%scatter plot of relative concentration and productions rate vs growth rate
figure(3)
plot(allFluorRel,allGrowthDeltaFluo,'ro','Markersize',5);
title('ASC825');
ylabel('Relative Growth Rate (dbl/hr)');
xlabel('Relative Production Rate (a.u.)');
ylim([-2 2]);
% xlim([-3000 3000]);
grid on;
hold on;
%%
%scatter plot conc vs growth rate
figure(4)
plot(allConc,allGrowthFluo,'bo','Markersize',1);
title('ASC1058');
ylabel('Relative Growth Rate (dbl/hr)');
xlabel('Relative Concentration (a.u.)');
ylim([-0.6 1]);
% xlim([100 350]);
grid on;
hold on;


