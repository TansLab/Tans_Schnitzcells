
%% 
%
% Preliminary analysis for switching experiments.

% Adapted code from MW_summaryplotspreliminaryanalysis.m
SCRIPTNAME = 'MW_summaryplotspreliminaryanalysisforswitch.m';

%% Parameters
% Dataset
if ~exist('MYDIR', 'var')
    %MYDIR = 'F:\A_Tans1_step1_incoming_not_backed_up\2016-02-10\outputSummary\';    
    %MYDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-05-15\outputSummary\';
    error('MYDIR var not set.')
end
if ~exist('SWITCHTIME', 'var')
    %SWITCHTIME = 171;
    %SWITCHTIME = (24*60+46)-(22*60+9)
    error('SWITCHTIME var not set.')
end
if ~exist('CUSTOMCOLORS', 'var')
    %CUSTOMCOLORS = linspecer(6);
    %CUSTOMCOLORS = linspecer(3);
    %CUSTOMCOLORS = linspecer(2);
    error('CUSTOMCOLORS var not set.')
end
if ~exist('FITWINDOWS', 'var')
    % give multiple fitwindows [window1t1, window1t2; window2..; etc]
    %FITWINDOWS = [0,SWITCHTIME; SWITCHTIME,800]; 
    error('FITWINDOWS var not set.')
end
if ~exist('FLUORIDXTOPLOT','var')
    error('Please set FLUORIDXTOPLOT');
    %FLUORIDXTOPLOT = [1]; % use numbers
    %FLUORIDXTOPLOT = [1 2]; % use numbers
end


% optional input
% MANUALSCHNITZLOCATIONS = {'..','..'};


%% Use MW_summaryplotspreliminaryanalysis to load data
if ~exist('flagMW_summaryplotspreliminaryanalysis','var')
    MW_summaryplotspreliminaryanalysis
end

%% Go over positions again and now create length/growth plots==============
h3=figure(); clf; 

subplot(1,2,1); hold on;

% Determine # positions
totalNrPositions = max(positionIndices);

% determine # fit windows
sizeFITWINDOWS = size(FITWINDOWS);
NrFitWindows = sizeFITWINDOWS(2);

% Plot summed length data--------------------------------------------------
fittedMus = nan(1,totalNrPositions); fittedL0=[];
for i = positionIndices
        
    % the data for all times
    time = [alldata(i).frameData(:).frameTime]./60; % convert to hours
    datatofit = log([alldata(i).frameData(:).frameSummedLength])./log(2); % fit exponential w. base 2
    
     % plot it
    l = plot([alldata(i).frameData(:).frameTime],[alldata(i).frameData(:).frameSummedLength],'o');
    set(l, 'LineWidth', 2, 'Color', CUSTOMCOLORS(i,:),'MarkerSize',10);    
    
    % fit it
    for j = 1:NrFitWindows
        % determine which data to fit (i.e. within fitting window)
        fitWindowHrs = FITWINDOWS./60;
        indexesInTimeWindow = time>fitWindowHrs(j,1) & time<fitWindowHrs(j,2);
        % fit it
        myfit = polyfit(time(indexesInTimeWindow),datatofit(indexesInTimeWindow),1);
        fittedMus(j,i) = myfit(1);
        fittedL0(j,i)  = 2^(myfit(2));

        % Plot fits as lines
        plot(fitWindowHrs(j,:)*60,fittedL0(j,i).*2.^(fittedMus(j,i)*fitWindowHrs(j,:)),...
            '-k','LineWidth',3)
    end
end

set(gca,'yscale','log');



% Determine ylim
frameSummedLengthdatapile=[];
for i=positionIndices
    frameSummedLengthdatapile = [frameSummedLengthdatapile, alldata(i).frameData(:).frameSummedLength];
end
theYLim = [min(frameSummedLengthdatapile)/2,max(frameSummedLengthdatapile)*2];

% Plot switchtime
plot([SWITCHTIME SWITCHTIME],theYLim,':k')

% Plot timewindow used for fittingfor j = 1:NrFitWindows
windowColors = linspecer(NrFitWindows);
for j = 1:NrFitWindows
    plot(FITWINDOWS(j,:),[theYLim(1) theYLim(1)],'-^','LineWidth',2,'MarkerFaceColor',...
            windowColors(j,:),'Color',windowColors(j,:));
end
%plot(FITWINDOWS(1),theYLim(1),'k>','LineWidth',2,'MarkerFaceColor','k')
%plot(FITWINDOWS(2),theYLim(1),'k<','LineWidth',2,'MarkerFaceColor','k')

% Set title etc.
title([SCRIPTNAME 10 'SWITCHTIME = ' num2str(SWITCHTIME) ' mins'],'Interpreter','None')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('summed length (a.u.)');
ylim(theYLim);

% bar plot of fitted growth rates------------------------------------------
subplot(1,2,2); hold on;

% ugly solution to issue that matlab treats n*1 matrix differently than n*m
if totalNrPositions==1
    fittedMus = [fittedMus, nan(NrFitWindows,1)];
end

% plot
b=bar(fittedMus');

% Set bar colors
for j = 1:NrFitWindows   
    set(b(j),'FaceColor', windowColors(j,:));
end


% plot
%{
for i=positionIndices
    l=bar(ones(1,numel(fittedMus(:,i)))*i, fittedMus(:,i));
    set(l, 'LineWidth', 1, 'FaceColor', CUSTOMCOLORS(i,:), 'EdgeColor', CUSTOMCOLORS(i,:));
end
%}

% align,labeling
set(gca, 'XTickLabel',{1:max(positionIndices)}, 'XTick',1:max(positionIndices))
xlim([0,max(positionIndices)+1])

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('colony number');
ylabel('fitted growth rate (doublings/hour)');

% save the figure 
saveas(h3, [MYDIR 'summaryPlot_colonylength_twowindows.tif'],'tif');
saveas(h3, [MYDIR 'summaryPlot_colonylength_twowindows.fig'],'fig'); %epsc

%% Summary plot with average of fitwindow1 mus vs. fitwindow2 mus

h4=figure(); clf; hold on;

% calculate means
notnanindices=~isnan(fittedMus(1,:));
windowmean(1) = mean(fittedMus(1,notnanindices));
windowmean(2) = mean(fittedMus(2,notnanindices));
windowstd(1) = std(fittedMus(1,notnanindices));
windowstd(2) = std(fittedMus(2,notnanindices));

% plotting
%[h,hErr] = barwitherr(windowstd,windowmean)
hErr = errorbar(windowmean,windowstd,'o')

set(hErr, 'LineWidth',3);

xlim([0,NrFitWindows+1]);
ylim([0,1.1]);

set(gca, 'XTickLabel',{1:NrFitWindows}, 'XTick',1:NrFitWindows)

xlabel('Fit window')
ylabel('Average growth rate +/- std [dbl/hr]')

MW_makeplotlookbetter(15);

saveas(h4, [MYDIR 'summaryPlot_colonylength_growthmean.tif'],'tif');
saveas(h4, [MYDIR 'summaryPlot_colonylength_growthmean.fig'],'fig');

%% Add switchtime to figure 1 (made in MW_summaryplotspreliminaryanalysis)

figure(h1); hold on;

plotrows=numel(FLUORIDXTOPLOT);

for fluorIdx = 1:plotrows

    subplot(plotrows,2,1+2*(fluorIdx-1)); hold on;
    % Plot switchtime
    plot([SWITCHTIME SWITCHTIME],fluorYlim1{fluorIdx},':k','LineWidth',2)

    subplot(plotrows,2,2+2*(fluorIdx-1)); hold on;
    % Plot switchtime
    plot([SWITCHTIME SWITCHTIME],fluorYlim2,':k','LineWidth',2)

end

saveas(h1, [MYDIR 'summaryPlot_colonylength_fluorsignal.tif'],'tif');
saveas(h1, [MYDIR 'summaryPlot_colonylength_fluorsignal.fig'],'fig');
