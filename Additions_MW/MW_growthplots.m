%% 2015/11 MW 
% 
% Plotting some growth data (DJK also made scripts for this..)


% Paramater settings ----------
if ~exist('FITTIME', 'var')
    FITTIME = [0 800];
end
if ~exist('EXPORTFOLDER', 'var')
    EXPORTFOLDER = 'C:\Users\wehrens\Desktop\export\';
end
if ~exist('FIGUREVISIBLE','var')
   FIGUREVISIBLE = 'on'; 
end
% -----------------------------


%%

timeField = 'time';
lengthField = 'length_fitNew';

%%

%%% CURRENT ONE: %%%
%myFile = '\\biofysica-store\data\TansGroup\Aileen\KILLERMIKE\2007mondate\2007-08-23\early Pos08-mini-01\data\Pos08-mini-01-Schnitz.mat';
if ~exist('myFile')
    myFile = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-06-02\pos4crop\data\pos4crop-Schnitz.mat';
end
%%% //////////// %%%

% create identifier
%slashes = find(myFile=='\');
%identifier = myFile(slashes(end-4):end);
%myID = regexprep(myFile,'\\','-');
%myID = regexprep(myID,'\.','-');

load(myFile); 
myFile % just print for user convenience


%% Analyze

myframe_nrs = unique([schnitzcells.frame_nrs]);
myTimes = unique([schnitzcells.(timeField)]);
allLengths = [schnitzcells.(lengthField)];

% Surely this can be done prettier, but let's just loop to get the desired
% data now.
lengthsPerframe_nrs = {}; plottingfXaxis = {}; timesPerframe_nrs = {};
for f = myframe_nrs
   
    lengthForThisFrame = [];
    timesForThisFrame = [];
    
    for schnitzIdx = 1:numel(schnitzcells)
        % Look whether this schnitz lives in current frame number,
        % and what timepoint belongs to the current frame.
        pointInSchnitz = find(schnitzcells(schnitzIdx).frame_nrs==f);
        
        % If so, add length at that timepoint to the collection of lengths
        if ~isempty(pointInSchnitz)
            timesForThisFrame(end+1) = schnitzcells(schnitzIdx).(timeField)(pointInSchnitz);            
            lengthForThisFrame(end+1) = schnitzcells(schnitzIdx).(lengthField)(pointInSchnitz);
        end
    end
    
    timesPerframe_nrs{end+1} = timesForThisFrame;
    lengthsPerframe_nrs{end+1} = lengthForThisFrame;
    plottingfXaxis{end+1} = ones(numel(lengthForThisFrame),1)*f; % convenient as x-axis later
    
end

%% Plotting

% Some plotting settings
h = figure('Name', myFile,'pos', [100 100 800+100 400+100],'Visible',FIGUREVISIBLE); % 700x600
subplot(1,2,1);
title([ 'Raw data, lengths bacteria at each frame'])
xlabel('Time in minutes');
ylabel('Length individual bacteria ({\mu}m)');
hold on;
xlim([min(myTimes),max(myTimes)])
ylim([min(allLengths),max(allLengths)])

summedLengthsPerFrame=[];
for i = 1:numel(plottingfXaxis)
    
    % raw data
    xdata = cell2mat(timesPerframe_nrs(i));
    ydata = cell2mat(lengthsPerframe_nrs(i));    
    plot(xdata, ydata,'o')
    set(gca,'FontSize',15)
        
    % determine sum of cell lengths
    summedLengthsPerFrame(end+1) = sum(cell2mat(lengthsPerframe_nrs(i)));
    
end

% Plotting setup
subplot(1,2,2); 
semilogy(myTimes, summedLengthsPerFrame,'o');
hold on;
ylim([min(summedLengthsPerFrame)*.5,max(summedLengthsPerFrame)*2]);

% Fitting
fit_idx = find(myTimes<FITTIME(2));
[fitMu A0] = DJK_ExponentialFit(myTimes(fit_idx)/60, summedLengthsPerFrame(fit_idx));

% Plotting
semilogy(myTimes,A0*2.^(fitMu*myTimes/60),'-k','LineWidth',2)
set(gca,'FontSize',15)

title([ 'Growth \newline' ...
        'Fitted \mu = ' num2str(fitMu) ' (dbl/hr)'])
xlabel('Time in minutes');
ylabel('Summed length all bacteria ({\mu}m)');

set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')

saveas(h, [EXPORTFOLDER 'EPS_bacteriaSize_RAW-' p.myID '.eps'],'epsc')
saveas(h, [EXPORTFOLDER 'TIF_bacteriaSize_RAW-' p.myID '.tif'])

%% Fancy plot for suppl. mat. article
% Plotting setup
FONTSIZE = 24;

h = figure('pos', [100 100 700+100 600+100],'Visible',FIGUREVISIBLE); % 700x600
semilogy(myTimes, summedLengthsPerFrame,'o','LineWidth',4,'Color',[.6 .6 .6], 'MarkerSize',12);
hold on;
ylim([min(summedLengthsPerFrame)/2,max(summedLengthsPerFrame)*2]);

% Fitting
fit_idx = find(myTimes<FITTIME(2));
[fitMu A0] = DJK_ExponentialFit(myTimes(fit_idx)/60, summedLengthsPerFrame(fit_idx));

%title([ 'Colony growth'])
xlabel('Time (min)');
ylabel('Log summed lengths for colony (a.u.)');
set(gca,'YTickLabel',[]);
set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal')
set(gca,'FontSize',FONTSIZE)

saveas(h, [EXPORTFOLDER 'EPS_bacteriaSize-' p.myID '.eps'],'epsc')
saveas(h, [EXPORTFOLDER 'TIF_bacteriaSize-' p.myID '.tif'])

% Plotting of fitted line:
semilogy(myTimes,A0*2.^(fitMu*myTimes/60),'--k','LineWidth',5)

saveas(h, [EXPORTFOLDER 'EPS_bacteriaSize_withFit-' p.myID '.eps'],'epsc')
saveas(h, [EXPORTFOLDER 'TIF_bacteriaSize_withFit-' p.myID '.tif'])

disp('Done.');



%%

%{
%% Using Schnitzcells code

% copied from Excel file
p1 = DJK_initschnitz('Pos09-mini-01crop','2008-02-23','e.coli.amolf','rootDir','F:\X_Other_datasets\2008mondate\', 'cropLeftTop',[1 1], 'cropRightBottom',[800 800],'fluor1','none','fluor2','none','fluor3','none','setup','setup1','softwarePackage','metamorph','camera','?');
FITTIME(2) = DJK_analyzeMu(p1, schnitzcells, 'xlim', [0 2000], 'onScreen', 0,'FITTIME(2)',[1000 1200]);
%}









