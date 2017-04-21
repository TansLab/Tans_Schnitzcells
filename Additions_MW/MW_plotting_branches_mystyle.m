function [h1, gatheredBranchOutput] = MW_plotting_branches_mystyle(p,schnitzcells,avgLineColor,whichParamsToPlot,fluorColorToPlot,theYLabels)
% [h1, gatheredBranchOutput] = MW_plotting_branches_mystyle(p,schnitzcells,avgLineColor,whichParamsToPlot,fluorColorToPlot,theYLabels)
%
% For this fn input
% - p                   p should have
%       p.fitTime
% - schnitzcells        standard schnitzcells struct
% - avgLineColor        color for avg line
% - whichParamsToPlot   {{'set1x','set1y'},..}
%                       e.g. {{'X_time',    'muP9_fitNew_cycCor'}},...
%                             {'dX5_time',  'dX5_cycCor'},...
%                             {'X_time',    'X6_mean_cycCor'}};
% - fluorColorToPlot    Y,C,G,R
% - theYLabels 
%                       e.g. {'Growth [dbls/hr]', 'Production fluor [a.u.]', 'Fluor concentration [a.u.]',}
%         
% output parameters
% - h1                      handle to figure
% - gatheredBranchOutput    values of e.g. mean line are given by: 
%                           gatheredBranchOutput(setIdx).binCenters,gatheredBranchOutput(setIdx).meanValuesForBins
%                           where setIdx corresponds with indices for whichParamsToPlot

%  User parameters
FONTSIZE=20;

%%

if ~exist('avgLineColor','var')
    avgLineColor='r';
end

if ~isfield(p,'fitTime')
    p.fitTime = [0,1000];
end

if ~exist('whichParamsToPlot','var')
    whichParamsToPlot = ...
                {{'X_time',    'muP9_fitNew_cycCor'},...
                 {'dX5_time',  'dX5_cycCor'},...
                 {'X_time',    'X6_mean_cycCor'}};
end

if ~exist('theYLabels','var')
    theYLabels = {'Growth [dbls/hr]',...
             'Production fluor [a.u.]'...
             'Fluor concentration [a.u.]',};
end

%%
        
for ii=1:numel(whichParamsToPlot)
    whichParamsToPlot{ii} = strrep(whichParamsToPlot{ii},'X',fluorColorToPlot);
    theYLabels{ii} = strrep(theYLabels{ii},'X',fluorColorToPlot);
end




%%
% note that the MW_delayedScatter script also has sophisticated branch
% plotting script.

%% gather data
gatheredBranchOutput=struct;
for setIdx=1:numel(whichParamsToPlot)

    %% Get branch data

    myDataFields = whichParamsToPlot{setIdx};

    p.dataFields = myDataFields;
    gatheredBranchOutput(setIdx).myBranches = DJK_get_branches(p,schnitzcells)
    
       
    % determine average line
    signalData={gatheredBranchOutput(setIdx).myBranches(:).(whichParamsToPlot{setIdx}{2})};
    timeData={gatheredBranchOutput(setIdx).myBranches(:).(whichParamsToPlot{setIdx}{1})};
    
    allTimes = unique([timeData{:}]);
    deltaTime = timeData{1}(2)-timeData{1}(1);
    myBins = [allTimes allTimes(end)+deltaTime]-deltaTime/2;
    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins] = binnedaveraging(timeData,signalData,myBins);
    
    gatheredBranchOutput(setIdx).binCenters  = binCenters;
    gatheredBranchOutput(setIdx).meanValuesForBins = meanValuesForBins;
end

%% plot it
HIGHLIGHTNR = 2;
MYSTYLES   = {'-',':',':'};

h1=figure(); clf; hold on; 

PLOTWIDTH=330;
set(h1, 'Position', [250 350 250+PLOTWIDTH*numel(whichParamsToPlot)  450]);

NrOfSubPlots = numel(whichParamsToPlot);

for setIdx=1:NrOfSubPlots
    
    %% Plot it  
    
    %subplot(1,numel(whichParamsToPlot),setIdx); hold on;
        
    %subtightplot(1,numel(whichParamsToPlot),setIdx,[0.08,0.08],0.15,0.07); hold on
    subtightplot(1,numel(whichParamsToPlot),setIdx,[0.08,0.24/NrOfSubPlots*1.5],0.01*FONTSIZE,0.012*FONTSIZE/NrOfSubPlots); hold on
    %subplot(1,numel(whichParamsToPlot),setIdx); hold on
    
    % plot branches
    totalNrBranches = numel(gatheredBranchOutput(setIdx).myBranches);
    for branchIdx = 1:totalNrBranches
        plot(gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{1}), gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{2}),'-','Color',[.6 .6 .6])
    end
    % color some randomly selected in different color
    colorBranchLinesIdx = ceil(rand(1,HIGHLIGHTNR)*(totalNrBranches));
    for ii = 1:numel(colorBranchLinesIdx)
        branchIdx = colorBranchLinesIdx(ii);
        %plot(gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{1}), gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{2}),'-','Color',[0 0 0],'LineWidth',2)
        plot(gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{1}), gatheredBranchOutput(setIdx).myBranches(branchIdx).(whichParamsToPlot{setIdx}{2}),...
            MYSTYLES{ii},'Color',[0 0 0],'LineWidth',1,'MarkerFaceColor','k','LineWidth',3)
    end
    

    % plot average line
    plot(gatheredBranchOutput(setIdx).binCenters,gatheredBranchOutput(setIdx).meanValuesForBins,'-','LineWidth',5,'Color',avgLineColor);
    
    % ylim
    myYlim=[0,max([gatheredBranchOutput(setIdx).myBranches(:).(whichParamsToPlot{setIdx}{2})])*1.1];
    ylim(myYlim);
    
    % xlim
    myXlim=[0, max(allTimes)];
    xlim(myXlim)
    
    % plot timeslot    
    plot([p.fitTime(1),p.fitTime(1)],myYlim,'-k');
    plot([p.fitTime(2),p.fitTime(2)],myYlim,'-k');
    plot([p.fitTime(1)],myYlim(1),'-^k','MarkerFaceColor','k');
    plot([p.fitTime(2)],myYlim(1),'-^k','MarkerFaceColor','k');
    
    % cosmetics   
    MW_makeplotlookbetter(FONTSIZE);
    xlabel('time [min]');
    ylabel(theYLabels{setIdx});        

end

end


















