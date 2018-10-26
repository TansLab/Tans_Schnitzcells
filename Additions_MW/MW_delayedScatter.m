
%% Description

% MW 2015/04
%
% NEW DESCRIPTION
% This script was historically used to generate delayed scatter plots, but
% has become a "key" script that performs cross-correlation analysis and
% also outputs other parameter information (like probability distrutions,
% branch plots, noise over time, etc). 
%     The script is not optimized for this task, and given the original
% function of calculating CCs, there is a high level of redundancy, which
% is especially inconvenient since the introduction of the controls.
%     For the CC analysis, but also for the others, the script starts with 
% the "p" and "schnitzcells" parameters, plus some additional parameters,
% to get CCs (and more). It is generally called from the 
% "Schnitzcells_masterscript.m" script.
% 
%
% OLD DESCRIPTION
% Script to quickly generate "delayed scatter" plots, for now based on
% NW's data handling; 
% MW TODO: Check PN data.
%
% About input params:
% - associatedFieldNames{1} MUST BE timefield field name
% - associatedFieldNames{2} X, usually fluor concentration or rate field name
% - associatedFieldNames{3} Y, usually growth rate field name
%
% - makeDtime   if set to 0 nothing happens, if set to 1, extra 
%               MWDJK_dX_time filed is recalculated.
% - myFitTime   This is the time that is taken into consideration for the
%               analysis. In Schnitzcells_masterscript, 
%               myFitTime = ourSettings.fitTimeCrosscorr;
%
% Additional (optional) parameters:
% - NOTMEANSUBTRACTED   : Don't automatically take fields that have the
%                         noise subtracted, i.e. noise_FieldOfInterest.
% - FIGUREVISIBLE       : choose 'off' or 'on' to make figures
%                         invisible/visible
% - p.recalcNoise       : re-calculates the colony means per frame and
%                         subtracts those
% - SHOWSOMECONTROLLINES  if this parameter is set then control line
%                          examples will be shown.
%
%
% TODO
% Note that it is actually a bit awkward that delayed scatter also
% performs analyses on the single parameter level (like plotting branches
% and plotting branch groups), since this leads to a load of redundancy.
% Similarly, a load of computation could be saved if branch groups are
% calculated for multiple parameters (>2) at once, such that multiple
% cross-correlations (or other two way comparisons) could be determined
% from one branch data set...
%
%
%
% Example of how to call script:
%{
CONFIGFILE = 'config_projectCRPcAMP'; % Not necessary this script, for later analysis scripts

myID = 'WT_pl-pRCRP-GFP_pl-CRP'; 
p.movieName = 'pos4crop';       % Not necessary if p already exists.
p.movieDate = '2015-06-12';     % Not necessary if p already exists.
p.fluor1='g';
myFitTime = [0 800]; 

associatedFieldNames =  {'G_time','G6_mean_cycCor', 'muP5_fitNew_cycCor'} % NW suggested fields

myIllumTime = 100; % Not necessary this script, for later analysis scripts
ASCnumber = 852; % Not necessary this script, for later analysis scripts
filterset='engfp'; % Not necessary this script, for later analysis scripts
fluoName = 'GFP'; % Not necessary this script, for later analysis scripts

ourSettings.myOutputFolder = ['F:\A_Tans1_step1_incoming_not_backed_up\'  p. movieDate   '\' p. movieDate  '_' p.movieName '_' myID  '\'];

p.NW_saveDir = [ourSettings.myOutputFolder 'misc\'];  % To send additional output to
p.DJK_saveDir = [ourSettings.myOutputFolder 'misc\']; % To send additional output to

% Location of .mat file containing schnitzcells struct
myDataFile = ['F:\A_Tans1_step1_incoming_not_backed_up\' p.movieDate '\' p.movieName  '\data\' p.movieName '-Schnitz.mat'];
load(myDataFile);

myTitle = 'WT CRP behavior July2, r1'; % Plot title

% info required to make branches
badSchnitzes = [868, 853, 774, 578]; % pos 4 bad ones

% Options for script
addSlowOnes = 1;            % Automatically add slow schnitzes to bad schnitzes
alreadyRemovedInMatFile=0;  % If s_rm is your input, it's not necessary to generate it, then set this to 1
makeDtime = 0;              % !! If 1, re-calculates schnitzcells.dX_time field !!
PLOTSCATTER=1;              % If zero, only cross corrs are calculated, not delayed scatter plots

% Run script
MW_delayedScatter.m 

% One can just set associatedFieldNames to another value, and re-run the
% script to obtain cross corrs etc. between other paramters.

%}

% some parameters for visuals
someColors=linspecer(5);
PLOTSIZE=[6,4]; FONTSIZE=7;

% ??
PERFORMSOMECHECKS = 0;

% additional options
if ~exist('myOutputFolder')
    myOutputFolder = 'C:\Users\wehrens\Desktop\testdelayedscatter\output\';
end

if ~exist('associatedFieldNames') | ~exist('p') | ~exist('badSchnitzes')
    error('input not supplied.')
end
    
if ~isfield(p,'sameLength')
    p.sameLength=1; % assume that the branches are the same length unless told otherwise
end

if ~isfield(p,'extraNorm')
    p.extraNorm=0;
end

if ~isfield(p,'recalcNoise')
    p.recalcNoise=1;
end


% At sections there are some more parameters to set. They are marked in
% capitals.

% Loading
%load(myDataFile); % Can also be done by user

%% 
% Collect indices of the GFP measurement per schnitz.
% (Second method)
% ===    
% This is rather specific conversion code, only use if you know what it is
% about; otherwise ignore. -MW 2015/08
% It grabs the values from timepoints corresponding to the timepoints where
% a fluor or dfluor was calculated for (this reference field is given by 
% FIELDOFREFERENCE).
%
% Fields that need to be set:
%{
FIELD1='time';
FIELD3='muP15_fitNew_all';
FIELDOFREFERENCE = 'G5_mean_all'; % also corresponds to dG5, except for branch end/start.
%}

if exist('RECALCULATE', 'var')
    if RECALCULATE == 1
        
        indicesForField = {};
        numelschnitzcells=numel(schnitzcells);
        for schIdx = 1:numelschnitzcells            
            
            % Get indices where myField has value
            indicesForField{schIdx} = find(~isnan(schnitzcells(schIdx).(FIELDOFREFERENCE))); % could access field name as string
            
            if ~(numel(schnitzcells(schIdx).(FIELDOFREFERENCE)) == numel(schnitzcells(schIdx).(FIELD1))) == ...
                   numel(schnitzcells(schIdx).(FIELD3))
                warning(['Things are messed up? Look at schIdx = ' num2str(schIdx)]);
            end

            schnitzcells(schIdx).([FIELD1 '_MW_atX']) = schnitzcells(schIdx).(FIELD1)(indicesForField{schIdx});
            schnitzcells(schIdx).([FIELD3 '_MW_atX']) = schnitzcells(schIdx).(FIELD3)(indicesForField{schIdx});
            
            % Now for dX
            % Assuming PN_fluorrate_X was used. (And FIELDOFREFERENCE
            % contains fluor value.)
            % First calculate times                
            dummySField = schnitzcells(schIdx).(FIELDOFREFERENCE);
            % if has no parent, remove value, like real value is also missing in PN_fluorrate_X
            if (schnitzcells(schIdx).P == 0) || ~any(~isnan(schnitzcells(schnitzcells(schIdx).P).(FIELDOFREFERENCE)))
                tempIdxs = find(~isnan(dummySField));
                if ~isempty(tempIdxs)
                    dummySField(tempIdxs(1))=NaN;
                end
            end
            % if has no daughter, remove value, like real value is also missing in PN_fluorrate_X
            if (schnitzcells(schIdx).E == 0) || ~any(~isnan(schnitzcells(schnitzcells(schIdx).E).(FIELDOFREFERENCE))) ...
                    || ~any(~isnan(schnitzcells(schnitzcells(schIdx).D).(FIELDOFREFERENCE)))
                tempIdxs = find(~isnan(dummySField));
                if ~isempty(tempIdxs)
                    dummySField(tempIdxs(end))=NaN;
                end
                %dummySField = dummySField(1:end-1);
            end
                        
            indicesForDField{schIdx} = find(~isnan(dummySField)); % could access field name as string
            schnitzcells(schIdx).([FIELD1 '_MW_atdX']) = schnitzcells(schIdx).(FIELD1)(indicesForDField{schIdx});
            schnitzcells(schIdx).([FIELD3 '_MW_atdX']) = schnitzcells(schIdx).(FIELD3)(indicesForDField{schIdx});
            
        end
    
    end
end

disp('done');

%% preparing data

name_rm = 'rm'; name_all = 'all';
fitTime = myFitTime;

if ~alreadyRemovedInMatFile
    
    % Find Schnitzes with slow/negative growth rate -> rm them?!
    slowschnitzes=NW_detectSlowSchnitzes(p,schnitzcells,associatedFieldNames{3},'muThreshold',0);

    % Now if badSchnitzes not given, just take the slow ones determined
    % above    
    if ~exist('badSchnitzes','var') 
        disp('ATTENTION: badSchnitzes not given, assuming slowsschnitzes are badSchnitzes..');
        badSchnitzes= slowschnitzes';
        pause(2);
    elseif exist('addSlowOnes','var'), if addSlowOnes
        
        disp('ATTENTION: Assuming slowsschnitzes are badSchnitzes, adding them to badSchnitzes array..');
        badSchnitzes=unique([badSchnitzes slowschnitzes']);
        pause(2);
    
    end, end
    
    % Preparation to load
    % Adapted from NW excel sheet    
    s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); 
    s_all_fitTime = DJK_selSchitzesToPlot(s_all, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_all_fitTime = ['all_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
    s_all_fitTime_cycle = DJK_selSchitzesToPlot(s_all_fitTime, 'completeCycle', @(x) x ~= 0); name_all_fitTime_cycle = [name_all_fitTime '_cycle'];

    s_rm = DJK_selSchitzesToPlot(s_all, 'P', @(x) 1); 
    if ~isempty(badSchnitzes)
        for branchIdx=badSchnitzes, s_rm(branchIdx).useForPlot=0; end;
    end
    s_rm_fitTime = DJK_selSchitzesToPlot(s_rm, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); name_rm_fitTime = ['rm_' num2str(fitTime(1)) '_' num2str(fitTime(2))];
    s_rm_fitTime_cycle = DJK_selSchitzesToPlot(s_rm_fitTime, 'completeCycle', @(x) x ~= 0); name_rm_fitTime_cycle = [name_rm_fitTime '_cycle'];
    warning('s_rm_fitTime_cycle is not used, this is a TODO, fix it! -MW');
end


%% Calculating branches
% ===
s_rm = MW_calculateframe_nrs(s_rm); % backwards compatibility fix 

fitTime = fitTime + [2 -2];

branchData = MW_getBranches(p,s_rm,'dataFields',{associatedFieldNames{1}, associatedFieldNames{2}, associatedFieldNames{3} }, 'fitTime', fitTime); 
    % set p.sameLength=0 to take all (unequally long) branches into account
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_Conc_oldRates'];


%% Plot branches
HIGHLIGHTSUSPICOUS = 0;
NBINS = 50;
%yfieldbranchtoplot = 3; % 3 = growth, 2 = fluor

for yfieldbranchtoplot=[2,3]

    %% Just some plot colors
    %distinguishableColors = distinguishable_colors(numel(branchData)+1,[1 1 1]);  % XXX
    % create distinguishable colors
    distinguishableColors = linspecer(numel(branchData)+1); 
    % randomize order of colors to get better effect
    distinguishableColors=distinguishableColors(randperm(size(distinguishableColors,1)),:);

    %% Plot all branches
    h1=figure('Visible', FIGUREVISIBLE); clf; hold on;
    offset=100; width1=800; height1=600;
    set(h1, 'Position', [offset offset width1 height1]);
    numelBranches = numel(branchData);
    branchLineHandles = [];
    branchDotsHandles = [];
    for branchIdx = 1:numelBranches
        branchLineHandles(end+1) = plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}),'-','Color',distinguishableColors(branchIdx,:));
        branchDotsHandles(end+1) = plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}),'.','Color',[1 1 1]);
        set(branchLineHandles(end), 'LineWidth', (numelBranches-branchIdx+1)/numelBranches*10);
    end

    %% Brute force average branches (note weighing is irrelevant for average)
    %{
    branchMatrixFieldX = []; branchMatrixFieldY = [];
    for branchIdx = 1:numelBranches
        currentXvector = branchData(branchIdx).(associatedFieldNames{1});
        branchMatrixFieldX = [branchMatrixFieldX; currentXvector];

        currentYvector = branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot});
        branchMatrixFieldY = [branchMatrixFieldY; currentYvector];
    end
    meanXvector = mean(branchMatrixFieldX);
    meanYvector = mean(branchMatrixFieldY);
    %}
    
    %% 
    currentXvector = {branchData(:).(associatedFieldNames{1})};
    currentYvector = {branchData(:).(associatedFieldNames{yfieldbranchtoplot})};
    uniqueTimes = unique([currentXvector{:}]);
    currentEdges = uniqueTimes-(uniqueTimes(2)-uniqueTimes(1))/2;
    [meanYvector, meanXvector,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianYvector]=binnedaveraging(currentXvector,currentYvector,currentEdges);    
    
    %%
    % xlabel
    xlabel(associatedFieldNames{1},'Interpreter', 'None'), ylabel(associatedFieldNames{yfieldbranchtoplot},'Interpreter', 'None')

    myXlimFig1 = max([branchData(:).(associatedFieldNames{1})]);
    xlim([0, myXlimFig1]);
    myYlimFig1 = [min([0, [branchData.(associatedFieldNames{yfieldbranchtoplot})]]),...
                  max([branchData.(associatedFieldNames{yfieldbranchtoplot})])];
%     ylim([myYlimFig1(1), myYlimFig1(2)*1.5]);


    %Set all fontsizes
    MW_makeplotlookbetter(20);
    
    %% Plot histogram
    h2=figure('Visible', FIGUREVISIBLE), clf, hold on
    allYdata = [branchData.(associatedFieldNames{yfieldbranchtoplot})];
    [nelements, centers] = hist(allYdata,NBINS)
    deltaY = centers(2)-centers(1);
    totalCount = numel(allYdata);
    %nelements=nelements./deltaY;
    %plot(centers,nelements,'or','LineWidth',2)
    pdfBarHandles = bar(centers,nelements,'FaceColor',someColors(1,:),'EdgeColor',someColors(1,:))    
    
    % Fit Normal distribution
    pd=fitdist(allYdata', 'Normal')    
    fittedDistrX = [min(allYdata):(max(allYdata)-min(allYdata))/100:max(allYdata)]
    fittedDistrYnorm = normpdf(fittedDistrX,pd.mu,pd.sigma)
    fittedDistrY = fittedDistrYnorm.*totalCount.*deltaY;
    % plot
    pdfLineHandle=plot(fittedDistrX,fittedDistrY,'k', 'LineWidth', 2);
    
    % Fit gamma distribution
    %{
    % Note that 
    % Hashimoto M, Nozoe T, Nakaoka H, Okura R, Akiyoshi S, Kaneko K, Kussell E, Wakamoto Y: Noise-driven growth rate gain in clonal cellular populations. Proc. Natl. Acad. Sci. U. S. A. 2016, 113:3251–3256.
    % claim the generation time has a gamma distribution; but we'd have to
    % transform the dbls/hr data -- or the distribution accordingly --
    % to check this claim. (Or make a plot of the generation times directly
    % from the schnitzcells.
    % Anyways, in this context it is not so useful.
    pdgamma=fitdist(allYdata(allYdata>0)', 'Gamma')    
    fittedDistrGammaX = [min(allYdata):(max(allYdata)-min(allYdata))/100:max(allYdata)];
    fittedDistrGammaYnorm = pdf(pdgamma,fittedDistrGammaX);
    fittedDistrGammaY = fittedDistrGammaYnorm .*totalCount.*deltaY;    
    plot(fittedDistrGammaX,fittedDistrGammaY,':k', 'LineWidth', 2);
    %}
    
    %probplot(allYdata)
    MW_makeplotlookbetter(20);

    title(['PDF for ' associatedFieldNames{yfieldbranchtoplot}],'Interpreter','None');
    xlabel(associatedFieldNames{yfieldbranchtoplot},'Interpreter','None');
    ylabel('PDF(s) * N * \Deltas (counts)');    

    % Define some axes limits
    myYlim = max(nelements)*1.1;
    ylim([0,myYlim]);
    xlim([myYlimFig1(1), myYlimFig1(2)*1.1]);

    % TODO make 99% confidence and plot in previous figure.
    %confidence = paramci(pd,'Alpha',.01);
    pdfSigmaHandles=[];
    sigma2 = pd.mu + 4.*[-pd.sigma, pd.sigma];
    pdfSigmaHandles(1)=plot([sigma2(1),sigma2(1)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    pdfSigmaHandles(2)=plot([sigma2(2),sigma2(2)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    sigma5 = pd.mu + 5.*[-pd.sigma, pd.sigma];
    pdfSigmaHandles(3)=plot([sigma5(1),sigma5(1)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)
    pdfSigmaHandles(4)=plot([sigma5(2),sigma5(2)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)

    %% Now also plot confidence intervals in previous figure
    %figure(h1); set(gcf,'Visible', FIGUREVISIBLE); hold on;
    sigmaLineHandles=[];
    set(0, 'CurrentFigure', h1); hold on;
    sigmaLineHandles(1)=plot([0,myXlimFig1],[sigma2(1),sigma2(1)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    sigmaLineHandles(2)=plot([0,myXlimFig1],[sigma2(2),sigma2(2)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    sigmaLineHandles(3)=plot([0,myXlimFig1],[sigma5(1),sigma5(1)],':','Color',[.5 .5 .5],'LineWidth', 2)
    sigmaLineHandles(4)=plot([0,myXlimFig1],[sigma5(2),sigma5(2)],':','Color',[.5 .5 .5],'LineWidth', 2)

    legend([sigmaLineHandles(1), sigmaLineHandles(3)],{'2\sigma confidence','5\sigma confidence'},'location','Best');

    %% Now list schnitzes that have suspiciously high signal:
    mySigma = sigma2; % i.e. the 2 sigma confidence interval
    suspiciousBranches = []; suspiciousSchnitzes = [];
    for branchIdx = 1:numel(branchData)
        if (     any(branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}) > mySigma(2)) ) || ...
           (     any(branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}) < mySigma(1)) )
            suspiciousBranches(end+1) = branchIdx;
            if HIGHLIGHTSUSPICOUS % plot if desired
                plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}),'-or','LineWidth',3)
            end
            %plot(branchData(branchIdx).(associatedFieldNames{1}),
            %branchData(branchIdx).(associatedFieldNames{YFIELDBRANCHPLOT}),'-o','LineWidth',3,'Color',mycolors(c)) % MW debug

            % Find out which schnitzes are suspiciously high or low
            locationsInThisBranch = unique(...
                [ find(branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}) > mySigma(2)),...
                  find(branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}) < mySigma(1)) ] ...
                );
            suspiciousSchnitzes = [suspiciousSchnitzes branchData(branchIdx).schnitzNrs(locationsInThisBranch)];
        end
    end
    suspiciousSchnitzes = unique(suspiciousSchnitzes)
    suspiciousBranches        

    %% Plot mean behavior
    %figure(h1); set(gcf,'Visible', FIGUREVISIBLE); hold on;
    set(0, 'CurrentFigure', h1); hold on;
    meanLineHandle=plot(meanXvector, meanYvector,'-','Color','k','LineWidth',3);
    medianLineHandle=plot(meanXvector, medianYvector,'-.','Color','k','LineWidth',3);
        
    %% Save the plots
    saveas(h1,[myOutputFolder 'TIF_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.tif']);
    saveas(h1,[myOutputFolder 'SVG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.svg']);
    saveas(h1,[myOutputFolder 'FIG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.fig']);
    %saveas(h1,[myOutputFolder 'EPS_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.eps'],'epsc');

    saveas(h2,[myOutputFolder 'TIF_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.tif']);
    saveas(h2,[myOutputFolder 'SVG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.svg']);
    saveas(h2,[myOutputFolder 'FIG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.fig']);
    %saveas(h2,[myOutputFolder 'EPS_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.eps'],'epsc');
    
    % Also save plots with better axes labels
    % figure(h1);  set(gcf,'Visible', FIGUREVISIBLE);
    set(0, 'CurrentFigure', h1); hold on;
    MW_makeplotlookbetter(FONTSIZE*2,[],PLOTSIZE);
    if associatedFieldNames{1,yfieldbranchtoplot}(1) == 'd'
        ylabel('Fluorophore production (a.u./min)');
    elseif associatedFieldNames{1,yfieldbranchtoplot}(1) == 'm'
        ylabel('Growth rate (doublings/hr)');
    elseif any(strcmp(upper(associatedFieldNames{1,yfieldbranchtoplot}(1)),{'G','R','C','Y'}))        
        ylabel('Fluorophore concentration (a.u./pixel)');
    end
    xlabel('Time (mins)');
    saveas(h1,[myOutputFolder 'TIF_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.tif']);
    saveas(h1,[myOutputFolder 'SVG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.svg']);
    saveas(h1,[myOutputFolder 'FIG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.fig']);
        
    %figure(h2); set(gcf,'Visible', FIGUREVISIBLE);
    set(0, 'CurrentFigure', h2); hold on;
    MW_makeplotlookbetter(FONTSIZE*2,[],PLOTSIZE);    
    if associatedFieldNames{1,yfieldbranchtoplot}(1) == 'd'
        xlabel('Production (a.u./min)');
    elseif associatedFieldNames{1,yfieldbranchtoplot}(1) == 'm'
        xlabel('Growth rate (doublings/hr)');
    elseif any(strcmp(upper(associatedFieldNames{1,yfieldbranchtoplot}(1)),{'G','R','C','Y'}))        
        xlabel('Concentration (a.u./pixel)');
    end
    title([]);
    ylabel('Observations');
    
    saveas(h2,[myOutputFolder 'TIF_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.tif']);
    saveas(h2,[myOutputFolder 'SVG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.svg']);
    saveas(h2,[myOutputFolder 'FIG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_readableLabels.fig']);
    
    %% Now modify the branch plot for smaller sized plots
    
    % Cosmetics
    set(0, 'CurrentFigure', h1); hold on;
    set(branchLineHandles, 'Color', [.8 .8 .8],'LineWidth',1)
    
    % Create highlight lines
    someMoreLineColors=linspecer(5);
    NRHIGHLIGHTLINES=3;
    for highlightLinesIdx=1:NRHIGHLIGHTLINES
        randomIdx=ceil(rand().*(numel(branchLineHandles)));
        set(branchLineHandles(randomIdx), 'Color', someMoreLineColors(highlightLinesIdx,:),'LineWidth',1);        
        uistack(branchLineHandles(randomIdx),'top');
    end
    
    % More cosmetics
    set(branchDotsHandles, 'Marker', 'none');
    set(meanLineHandle, 'LineWidth', 1);
    set(medianLineHandle, 'LineWidth', 1);
    set(sigmaLineHandles, 'Color', 'k','LineWidth',1);
    uistack(meanLineHandle,'top');
    legend('off');  
    
    % Y limits
    myYlimBranches = [min([branchData(:).(associatedFieldNames{yfieldbranchtoplot}),sigma5(1)]),...
              max([branchData(:).(associatedFieldNames{yfieldbranchtoplot}),sigma5(2)]) ];
    dataRange = myYlimBranches(2)-myYlimBranches(1);
    myYlimBranches=myYlimBranches+[-dataRange/100 dataRange/100]; %widen range by 1 percent
    ylim(myYlimBranches);
    
    % Front size
    MW_makeplotlookbetter(10);
    
    saveas(h1,[myOutputFolder 'TIF_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_small.tif']);
    saveas(h1,[myOutputFolder 'SVG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_small.svg']);
    saveas(h1,[myOutputFolder 'FIG_branches_' associatedFieldNames{1,yfieldbranchtoplot} '_small.fig']);
    
    %% Also for PDF create smaller version
    set(0, 'CurrentFigure', h2); hold on;
    if ishandle(pdfBarHandles), delete(pdfBarHandles); end
      
    pdfAreaHandles = area(centers,nelements,'EdgeColor','none');
    set(pdfAreaHandles,'FaceColor',[.7 .7 .7]);
    uistack(pdfAreaHandles,'bottom');
        
    set(pdfSigmaHandles, 'Color', 'k','LineWidth',1);
    set(pdfLineHandle, 'Color', 'k','LineWidth',1);
    
    MW_makeplotlookbetter(10);
    
    saveas(h2,[myOutputFolder 'TIF_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_small.tif']);
    saveas(h2,[myOutputFolder 'SVG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_small.svg']);
    saveas(h2,[myOutputFolder 'FIG_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '_small.fig']);
    %}
    %% For later output
    output.branchavg.([associatedFieldNames{yfieldbranchtoplot} '_xfield']) = meanXvector; % usually time
    output.branchavg.(associatedFieldNames{yfieldbranchtoplot}) = meanYvector;
    output.branchdata.(associatedFieldNames{yfieldbranchtoplot}) = branchData;
    %output.branchgroups.(associatedFieldNames{yfieldbranchtoplot}) = branch_groups;
    %output.branchgroups.(associatedFieldNames{yfieldbranchtoplot}) = branch_groupsControl;
    
    % save data relating to the pdf
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).centers           =centers;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).nelements         =nelements;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).fittedDistrX      =fittedDistrX;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).fittedDistrYnorm  =fittedDistrYnorm;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).fittedDistrY      =fittedDistrY;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).sigma2            =sigma2;
    output.pdf.(associatedFieldNames{yfieldbranchtoplot}).sigma5            =sigma5;
        
    
end


  
%% Get the actual cross-corrs

%REDUNDANCYALLOWED = 2^2;
REDUNDANCYALLOWED = 2^2;
if strcmp(FIGUREVISIBLE,'off')
    ONSCREEN=0;
else
    ONSCREEN=1;
end
NRBRANCHGROUPS=4;
FIELDPREFIX = 'noise_';
CONTROLSUFFIX = '_randomizedlineages';
%FIELDPREFIX = 'relative_';

% Some additional editing of the branches:
branchData = DJK_addToBranches_noise(p, branchData,'dataFields',{associatedFieldNames{1},associatedFieldNames{2},associatedFieldNames{3}});
if p.recalcNoise
    branchData = MW_addToBranches_noise(branchData, associatedFieldNames);
end

% Trim branches if they are the same length
if p.sameLength==1
    % Note that this only work if branches all originate at t=0, which
    % might not be the case if they are not the same length
    
    % Trim of starting frames until there are N start schnitzes
    trimmed_branches = DJK_trim_branch_data(branchData,NRBRANCHGROUPS);
    % Divide branchdata in groups based on those N start schnitzes
    branch_groups = DJK_divide_branch_data(trimmed_branches);
else    
    % Remove empty branches (can be case due short branches without points
    % where data was taken.
    branchData = MW_remove_empty_branches(branchData,3); 
    % Put in one group just to accomodate the script
    branch_groups = struct;
      
    branch_groups(1).branches = branchData;
    branch_groups(1).parent_cell = 0;
    branch_groups(1).nr_branches = length(branchData);
end

% In case you want to skip the branch group procedure, and create branch
% groups that are in fact just the original branches
% branchData = MW_remove_empty_branches(branchData,3);
% branch_groups = DJK_divide_branch_data(branchData);

% Colony average mean has already been substracted so theoretically extra
% normalization shouldn't have an effect.
%p.extraNorm=0;

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
[CorrData,composite_corr] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, [FIELDPREFIX associatedFieldNames{1,2}],[FIELDPREFIX associatedFieldNames{1,3}] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 

disp('Done getting cross correlations');

%% Create negative control 
clear multipleCorrDataControl
NUMCONTROLS=50;

% Create a control branch structure
if p.sameLength==1
    %% For negative control, combine growth rates and fluor signal traces randomly    
    % ===
    
    multipleCorrDataControl=struct; multipleComposite_corrControl=struct;
    if isfield(p,'dontmakeplots'), oldSettingdontmakeplots = p.dontmakeplots; else, oldSettingdontmakeplots=0; end
    p.dontmakeplots=1;
    % actually, get 100 negative controls here..    
    for repeatIdx= 1:NUMCONTROLS

        if mod(repeatIdx,10)==0
            disp(['Running controls, now at ' num2str(repeatIdx) '/' num2str(NUMCONTROLS) '.']);
        end
        
        % Suppress too many messages after 1st loop
        if repeatIdx==2
            warning('off')
        end
        
        % First initialize copy
        branch_groupsControl = branch_groups;

       % create a list like dataIdx(i) = [branchGroupIdx; branchIdx] that 
       % tells us where a certain branch i is located
       branchLocations=[];
       for grIdx = 1:numel(branch_groupsControl)
            branchLocations = [branchLocations ...            
                [grIdx.*ones(1,numel(branch_groupsControl(grIdx).branches));
                1:numel(branch_groupsControl(grIdx).branches)]
            ];
       end

       % now create a random permutation of these branches
       reshuffleIdxs = randperm(size(branchLocations,2));

       % Now reshuffle the branches
       for brIdx = 1:numel(reshuffleIdxs)
            randomSourceIdx = branchLocations(:,reshuffleIdxs(brIdx));
            branch_groupsControl(branchLocations(1,brIdx)).branches(branchLocations(2,brIdx)).([FIELDPREFIX associatedFieldNames{3} CONTROLSUFFIX]) = ...
                branch_groupsControl(randomSourceIdx(1)).branches(randomSourceIdx(2)).([FIELDPREFIX associatedFieldNames{3}]);
       end

        %{
        % For negative control, combine growth rates and fluor signal traces randomly
        % control one becomes equal to original
        branch_groupsControl = branch_groups;    
        for groupIdx=1:numel(branch_groups)

            % Clear last field (which is field that's going to be scrambled)
            branch_groupsControl(groupIdx).(associatedFieldNames{1,3}) = [];
            branch_groupsControl(groupIdx).([FIELDPREFIX associatedFieldNames{1,3}]) = [];

            % Get original data for last field
            data = {branch_groups(groupIdx).branches(:).(associatedFieldNames{1,3})};
            % Randomize data for last field
            randomizeddata = {data{randperm(numel(data))}};

            % Repeat above for noise_(..) fields
            % Get original data for last field
            noisedata = {branch_groups(groupIdx).branches(:).([FIELDPREFIX associatedFieldNames{1,3}])};
            % Randomize data for last field
            noiserandomizeddata = {noisedata{randperm(numel(noisedata))}};

            % for each branch
            for i=1:numel(branch_groups(groupIdx).branches)

                % set to randomized lineages within that branch
                branch_groupsControl(groupIdx).branches(i).([associatedFieldNames{1,3} CONTROLSUFFIX]) = ...
                    randomizeddata{i};

                % set to randomized lineages within that branch
                branch_groupsControl(groupIdx).branches(i).([FIELDPREFIX associatedFieldNames{1,3} CONTROLSUFFIX]) = ...
                    noiserandomizeddata{i};
            end
        end
        %}
        % Calculate the cross-correlation
        [CorrDataControl, composite_corrControl] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groupsControl, [FIELDPREFIX associatedFieldNames{1,2}],[FIELDPREFIX associatedFieldNames{1,3} CONTROLSUFFIX] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 
    
        % save them
        multipleCorrDataControl(repeatIdx).CorrDataControl                = CorrDataControl;
        multipleComposite_corrControl(repeatIdx).composite_corrControl    = composite_corrControl;
        
    end
    
    % re-activate warnings
    warning('on'); 
    % re-establish old plotting setting    
    p.dontmakeplots=oldSettingdontmakeplots;
else
    multipleCorrDataControl=struct; multipleComposite_corrControl=struct;
    if isfield(p,'dontmakeplots'), oldSettingdontmakeplots = p.dontmakeplots; else, oldSettingdontmakeplots=0; end
    p.dontmakeplots=1;
    % actually, get 100 negative controls here..
    for repeatIdx= 1:NUMCONTROLS

        if mod(repeatIdx,10)==0
            disp(['Running controls, now at ' num2str(repeatIdx) '/' num2str(NUMCONTROLS) '.']);
        end
        
        % Suppress too many messages after 1st loop
        if repeatIdx==2
            warning('off')
        end
        
        % when we have differently sized branches, the above control becomes
        % impossible to create -- so we resort to a less stringent control,
        % where at each point in time, we just randomly scramble the data from
        % all branches for one parameter.
        % Note that in this situation, there should only be 1 branchgroup, so
        % only branch_groupsControl(1) will exist. (See above why.)
        %%
        branch_groupsControl = branch_groups;
        %branch_groupsControl(1).branches(1).([associatedFieldNames{1,3} CONTROLSUFFIX]) = [];
        %branch_groupsControl(1).branches(1).([FIELDPREFIX associatedFieldNames{1,3} CONTROLSUFFIX]) = [];

        allTimePoints = unique([branch_groupsControl(1).branches.(associatedFieldNames{1})]);
        timeLookupTable={};

        for currentTimePointIdx=1:numel(allTimePoints)
            %%
            % get applicable indices for this point in time
            currentTimePoint=allTimePoints(currentTimePointIdx);

            % Collect indices for each branch that match with this timepoints
            % (empty if timepoint does not exist in this data)
            timeHits= arrayfun(@(br) find(branch_groupsControl(1).branches(br).(associatedFieldNames{1}) == currentTimePoint), 1:numel(branch_groupsControl(1).branches),'UniformOutput',0);
            % Branches for which there is data
            hitBranches= find(arrayfun(@(br) ~isempty(find(branch_groupsControl(1).branches(br).(associatedFieldNames{1}) == currentTimePoint)), 1:numel(branch_groupsControl(1).branches)));
            % Create lookup tables
            %timeLookupTable{currentTimePointIdx} = timeHits;
            %branchLookupTable{currentTimePointIdx} = hitBranches;

            % Now we know which branches to mix, so do that
            randomizedBranches = hitBranches(randperm(numel(hitBranches)));

            for brIdx=1:numel(hitBranches)

                sourceBranchIdx = hitBranches(brIdx);
                targetBranchIdx = randomizedBranches(brIdx);

                % actual randomization            
                branch_groupsControl(1).branches(targetBranchIdx).([associatedFieldNames{3} CONTROLSUFFIX])(timeHits{targetBranchIdx}) = ...
                    branch_groups(1).branches(sourceBranchIdx).([associatedFieldNames{3}])(timeHits{sourceBranchIdx});
                
                % also do for normalized field
                branch_groupsControl(1).branches(targetBranchIdx).([FIELDPREFIX associatedFieldNames{3} CONTROLSUFFIX])(timeHits{targetBranchIdx}) = ...
                    branch_groups(1).branches(sourceBranchIdx).([FIELDPREFIX associatedFieldNames{3}])(timeHits{sourceBranchIdx});

            end

        end

        % Calculate the cross-correlation
        [CorrDataControl, composite_corrControl] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groupsControl, [FIELDPREFIX associatedFieldNames{1,2}],[FIELDPREFIX associatedFieldNames{1,3} CONTROLSUFFIX] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 
        
        % save them
        multipleCorrDataControl(repeatIdx).CorrDataControl          = CorrDataControl;
        multipleComposite_corrControl(repeatIdx).composite_corrControl    = composite_corrControl;
        
    end
    % re-activate warnings
    warning('on'); 
    % re-establish old plotting setting    
    p.dontmakeplots=oldSettingdontmakeplots;
end

%% Now calculate the delayed scatter plots

% Do we want to filter out colony average behavior for the "delayed
% scatter" plots also? Maybe do this with noise fields?
% But let's try with "raw" data first..
p.timeField = associatedFieldNames{1,1};
% p.tauIndices = [-7:7]; %p.tauIndices = [-29:4:-1,0,1:4:30];
if isfield(p,'tauIndices'), p=rmfield(p,'tauIndices'); end
[dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = ...
    MW_getdelayedscatter(p, branchData, [FIELDPREFIX associatedFieldNames{1,2}], [FIELDPREFIX associatedFieldNames{1,3}], REDUNDANCYALLOWED);

% simple plot:
% halfwayidx = (numel(iTausCalculated)+1)/2; figure; scatter(dataPairsPerTau{halfwayidx}(:,1), dataPairsPerTau{halfwayidx}(:,2))

%% 
if ~exist('NOTMEANSUBTRACTED','var')
    NOTMEANSUBTRACTED = 0;
end
if NOTMEANSUBTRACTED
    [dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = ...
    MW_getdelayedscatter(p, branchData, [associatedFieldNames{1,2}], [associatedFieldNames{1,3}], REDUNDANCYALLOWED)
end

%% Plot trace of mean point..
if NOTMEANSUBTRACTED
    
    % Initialize
    if ~exist('SHOWPLUSMINFROMZERO','var'), SHOWPLUSMINFROMZERO = 25; end
    
    % Set up figure
    h98=figure('Visible', FIGUREVISIBLE); clf; hold on;
    
    % Redundant with below
    indexMidpoint = ceil(numel(iTausCalculated)/2)
    rangeiTausCalculated = [indexMidpoint-SHOWPLUSMINFROMZERO:indexMidpoint+SHOWPLUSMINFROMZERO];
    
    % Loop over different values of tau (= delays)
    myTimeTrace = [];
    numelRangeiTausCalculated = numel(rangeiTausCalculated);
    for groupIdx = 1:numelRangeiTausCalculated
        % Create array with the points.
        myTimeTrace = [myTimeTrace; mean(dataPairsPerTau{rangeiTausCalculated(groupIdx)}(:,1)) , mean(dataPairsPerTau{rangeiTausCalculated(groupIdx)}(:,2))];
        
        % Plot separately, color coded for amount of delay        
        l = plot(myTimeTrace(groupIdx,1), myTimeTrace(groupIdx,2),'o');
        % Color of plot
        timeColor = [0 0 1-groupIdx/numelRangeiTausCalculated]  + ... % Starting with blue
                    [groupIdx/numelRangeiTausCalculated 0 0];       % turning red, over time        
        set(l, 'Color', timeColor, 'MarkerFaceColor', timeColor);
            
    end
    
    % Plot time trace 
    plot(myTimeTrace(:,1), myTimeTrace(:,2),'-k')
    
    % 
    xlim([  min(dataPairsPerTau{indexMidpoint}(:,1)), max(dataPairsPerTau{indexMidpoint}(:,1))  ])
    ylim([  min(dataPairsPerTau{indexMidpoint}(:,2)), max(dataPairsPerTau{indexMidpoint}(:,2))  ])
    
    title(['< X(t) Y(t+\tau)>_{t}; \tau trace ' ...
        10 'Branches assumed independent, means shouldn''t move.'])    
        %10 'I''m not sure what this plot means; or should behave like??' ...
        %10 'this might all be an artefact of data selection..'])
    xlabel('<X(t)>_t');
    ylabel('<Y(t+\tau)>_t');
        
    MW_makeplotlookbetter(15);
end

%% Plot "raw" cross cor I calculate (MW)

myfig=figure('Visible', FIGUREVISIBLE); clf; hold on;
l=plot(iTausCalculated,correlationsPerTau,'o-r','LineWidth',2);

%% Compare two cross-corrs (DJK & MW), also plot the control
h5=figure('Visible', FIGUREVISIBLE); clf; hold on;

myColorsLS = linspecer(4); myColors = [0 0 0; myColorsLS(2,:); myColorsLS(1,:)];
%myColors = [0 0 0; 1 0 0; 0 70/255 170/255];

% Plot control DJK cross correlation function
% ===
%errorbar(CorrDataControl(:,1),CorrDataControl(:,2),CorrDataControl(:,3),'x-','Color', [.5,.5,.5], 'LineWidth',2)
    % myColorControl=[0 204 255]/255;  myColorControl=[0 85 212]/255;   myColorControl = myColors(3,:);
    myColorControl = myColors(3,:);
if exist('multipleCorrDataControl','var')
    
    % determine the boundaries of the area that is touched by the control
    % lines
    timeDataControl= multipleCorrDataControl(1).CorrDataControl(:,1)';
    myLineMeanControl=nan(1,numel(timeDataControl)); myLineMax=nan(1,numel(timeDataControl)); myLineMin=nan(1,numel(timeDataControl));
    for idx=1:numel(timeDataControl)
        
        % determine the line for each timepoint
        valuesAtTau=arrayfun(@(idx2) [multipleCorrDataControl(idx2).CorrDataControl(idx,2)], 1:numel(multipleCorrDataControl));
        myLineMeanControl(idx) = mean(valuesAtTau);
        myLineMax(idx)  = max(valuesAtTau);
        myLineMin(idx)  = min(valuesAtTau);
        
        %l3=plot(multipleCorrDataControl(idx).CorrDataControl(:,1),multipleCorrDataControl(idx).CorrDataControl(:,2),...
        %    '-','LineWidth',2,'Color',[.7 .7 .7])
    end
    
    % Now show some example lines if desired
    if exist('SHOWSOMECONTROLLINES','var')
        lighterShadeColor=[.7 .7 .7]; %lighterShadeColor=min(1,[((myColorControl+.7))]);
        for idx= 1:3%numel(multipleCorrDataControl)

            lSc=plot(multipleCorrDataControl(idx).CorrDataControl(:,1),multipleCorrDataControl(idx).CorrDataControl(:,2),...
                    '-','LineWidth',1,'Color', lighterShadeColor)

            % highlight some lines
            %{
            if any(idx==[numel(multipleCorrDataControl)-2:numel(multipleCorrDataControl)])
                set(lSc, 'Color', myColorControl);
            end
            %}
        end
    end
    
   %% Plot the control
   lControlArea1=area(timeDataControl,myLineMax,...
       'FaceColor',myColorControl,'LineStyle','none','FaceAlpha',.25);
   %alpha(fh,.25)
   lControlArea2=area(timeDataControl,myLineMin,...
       'FaceColor',myColorControl,'LineStyle','none','FaceAlpha',.25);
   %alpha(fh,.25)
   
   lControlMean=plot(timeDataControl,myLineMeanControl,'-',...
       'LineWidth',2,'Color',myColorControl);
      
   %{
   plot(timeData,myLineMax,'-',...
       'LineWidth',2,'Color',myColorControl);
   plot(timeData,myLineMin,'-',...
       'LineWidth',2,'Color',myColorControl);
    %}
    %{
    for idx=1:numel(multipleCorrDataControl)
        lControlMean=plot(multipleCorrDataControl(idx).CorrDataControl(:,1),multipleCorrDataControl(idx).CorrDataControl(:,2),...
            '-','LineWidth',2,'Color',[.7 .7 .7])
    end
    %}
else
    lControlMean=plot(CorrDataControl(:,1),CorrDataControl(:,2),'x-','LineWidth',2,'Color',myColorControl)
end

%% Plot MW scatter cross correlation function
% ===
% Calculate appropriate x-axis assuming assuming same delta(x) as CorrData,
% and dx is same everywhere.
centerIdx=ceil(size(CorrData,1)/2);
deltaXCorrData = CorrData(centerIdx+1,1)-CorrData(centerIdx,1);
MWxAxis = iTausCalculated.*deltaXCorrData;

% PLOT make the actual plot of scatter correlations
lCCmw=plot(MWxAxis,correlationsPerTau,'-','LineWidth',1,'MarkerFaceColor',myColors(2,:),'Color',myColors(2,:));
lCCdotsmw=plot(MWxAxis,correlationsPerTau,'.','LineWidth',1,'MarkerFaceColor',myColors(2,:),'Color',myColors(2,:));

%% Plot DJK cross correlation function
% ===
if size(CorrData,2)>2
    % error bars from branchgroups
    hCCerrorDJK=errorbar(CorrData(:,1),CorrData(:,2),CorrData(:,3),'-','Color', myColors(1,:), 'LineWidth',1)
    hCCerrorDJK.CapSize = 1;
end
% line itself
lCCdjk=plot(CorrData(:,1),CorrData(:,2),'-','LineWidth',2,'MarkerFaceColor',myColors(1,:),'Color',myColors(1,:))
lCCdotsDjk=plot(CorrData(:,1),CorrData(:,2),'.','LineWidth',1,'MarkerFaceColor',[.5 .5 .5],'Color',[.5 .5 .5])
%l1=plot(CorrData(:,1),CorrData(:,2),'.','LineWidth',2,'MarkerFaceColor',[1 1 1],'Color',[1 1 1])

% cosmetics
% ===

% If you recalculate correlations again w. different params, this allows
% plotting of extra line.
%l3=plot(CorrData(:,1),correlationsPerTau100,'o-','Color',[1 .5 0],'LineWidth',2)

myxlimvalues=[min(CorrData(:,1)), max(CorrData(:,1))];
xlim(myxlimvalues);
%ylim([-1.1,1.1]);
if ~strcmp(associatedFieldNames{1,2},associatedFieldNames{1,3})
    myylimvalues=[min(CorrData(:,2)), max(CorrData(:,2))];
    if myylimvalues(1)>-.5, myylimvalues(1) = -.5; end
    if myylimvalues(2)<.5, myylimvalues(2) = .5; end
    ylim(myylimvalues); % normal cross-correlation usually within this range

else
    myylimvalues=[min(CorrData(:,2)), 1];
    if myylimvalues(1)>-.2, myylimvalues(1) = -.2; end
    ylim(myylimvalues); % when fields are identical, it is an autocorrelation, and limits should be adjusted accordingly
end

plot(myxlimvalues,[0,0],'k-');

%legend([l1,l2,l3],{'DJK','MW','Control'})
%legend([l1,l2,l3],{'Lineages','Scatters','Control'})
legend([lCCdjk,lCCmw,lControlMean],{'Lineages','Scatters','Control'});%,'Location','Best')

if p.sameLength==0, sameLengthComment=[10 ' Composite of unequally-sized branches']; else sameLengthComment=''; end
title(['DJK vs. MW ' 10 myID '_' p. movieDate  '_' p.movieName sameLengthComment], 'Interpreter', 'none');
xlabel('\tau (hrs)');
ylabel(['R(' associatedFieldNames{1,2} ',... ' 10 associatedFieldNames{1,3} ') (normalized)'], 'Interpreter', 'none');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

plot([0,0],[-1,1],'-k');

saveas(h5,[myOutputFolder 'TIF_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '.tif']);
saveas(h5,[myOutputFolder 'EPS_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '.eps'],'epsc');

%% Now make it look nicer for saving

%figure(h5); set(gcf,'Visible', FIGUREVISIBLE);
set(0, 'CurrentFigure', h5); hold on;
MW_makeplotlookbetter(FONTSIZE*2,[],PLOTSIZE);    
title([]);
xlabel('Delay (hrs)');
legend([lCCdjk,lCCmw,lControlArea1],{'Lineages','Scatters','Control'},'Location','Best')
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

% Determine y label
% ===
% Simply based on the parameter name
if associatedFieldNames{1,2}(1) == 'd'
    param1name='p';
elseif associatedFieldNames{1,2}(1) == 'm'
    param1name='µ';
elseif any(strcmp(upper(associatedFieldNames{1,2}(1)),{'G','R','C','Y'}))        
    param1name='C';
end
% Again for 2nd parameter
if associatedFieldNames{1,3}(1) == 'd'
    param2name='p';
elseif associatedFieldNames{1,3}(1) == 'm'
    param2name='µ';
elseif any(strcmp(upper(associatedFieldNames{1,3}(1)),{'G','R','C','Y'}))        
    param2name='C';
end
ylabel(['Correlation(' param1name ',' param2name ')']);

%
saveas(h5,[myOutputFolder 'FIG_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_humanReadable.fig']);
saveas(h5,[myOutputFolder 'TIF_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_humanReadable.tif']);
saveas(h5,[myOutputFolder 'SVG_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_humanReadable.svg']);
%saveas(h5,[myOutputFolder 'EPS_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '.eps'],'epsc');

%% Now also modify it for smaller visualization
set(0, 'CurrentFigure', h5); hold on;

set(lControlMean, 'LineWidth',1,'Color',myColorControl);
set([lControlArea1 lControlArea2], 'FaceColor',[.8 .8 .8],'EdgeColor','none','FaceAlpha',1);
set(lCCdjk, 'LineWidth',1,'Color','k');
legend('off');

NRERRORBARS=20;
selectedIndicesError=round(linspace(1,numel(CorrData(:,1)),NRERRORBARS));
%delete(hCCerrorDJKsmall)
hCCerrorDJKsmall=errorbar(CorrData(selectedIndicesError,1),CorrData(selectedIndicesError,2),CorrData(selectedIndicesError,3),'LineStyle','none','Color','k','CapSize',0)

if ishandle(hCCerrorDJK), delete(hCCerrorDJK); end
if ishandle(lCCdotsDjk), delete(lCCdotsDjk); end
if ishandle(lCCdotsmw), delete(lCCdotsmw); end

MW_makeplotlookbetter(10);

saveas(h5,[myOutputFolder 'FIG_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_small.fig']);
saveas(h5,[myOutputFolder 'TIF_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_small.tif']);
saveas(h5,[myOutputFolder 'SVG_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_small.svg']);

%% Now also save these CC values for later usage

% CC based on scatter plots
output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).iTausCalculated = iTausCalculated;    
output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).correlationsPerTau = correlationsPerTau;
output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).MWxAxis = MWxAxis;

% Control CC
if exist('multipleCorrDataControl','var')

    output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).multipleCorrDataControl = multipleCorrDataControl;
    
    output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).myLineMeanControl = myLineMeanControl;
    output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).myLineMax = myLineMax;
    output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).myLineMin = myLineMin;

else
    
    output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).CorrDataControl = CorrDataControl;
    
end

% CC itself
output.CC.([associatedFieldNames{1,2} '_' associatedFieldNames{1,3}]).CorrData = CorrData;

%% Plot code from CRPcAMP..overview..general
% ==========
if ~exist('PLOTSCATTER','var'), PLOTSCATTER=1; end;
if PLOTSCATTER
NRCONTOURLINES = 5;
SHOWPLUSMINFROMZERO = 25;
PLOT3DSCATTER = 0;

% Whether plots should be shown.
if ~exist('FIGUREVISIBLE','var'), FIGUREVISIBLE='on'; end;

% What range should be plotted.
numeliTausCalculated=numel(iTausCalculated);
indexMidpoint = ceil(numeliTausCalculated/2)
ShowPlusMinFromZero=min(SHOWPLUSMINFROMZERO, indexMidpoint-1);
rangeiTausCalculated = [indexMidpoint-ShowPlusMinFromZero:indexMidpoint+ShowPlusMinFromZero];
numelRangeiTausCalculated = numel(rangeiTausCalculated);
% delayIdx = 11; % 39 is middle

myColorMap = colormap(winter(numel(iTausCalculated)));

if PLOT3DSCATTER
    h2 = figure('Visible', FIGUREVISIBLE); clf; hold on;
    offset=100; width1=500; height1=500;
    set(h2, 'Position', [offset offset width1 height1]);
end

h3 = figure('Visible', FIGUREVISIBLE); 
clf; hold on;

% Initialization
emptySizedCell = cell(1,numeliTausCalculated);
bandwidths = emptySizedCell; densities=emptySizedCell; 
Xs=emptySizedCell; Ys=emptySizedCell; Zs=emptySizedCell; count = 0;
for delayIdx = rangeiTausCalculated
    
    % Rename data more convenient
    data = [dataPairsPerTau{delayIdx}(:,1), dataPairsPerTau{delayIdx}(:,2)];    
    
    % plot scatter
    if PLOT3DSCATTER
        %h2 = figure(h2); set(gcf,'Visible', FIGUREVISIBLE);
        set(0, 'CurrentFigure', h2); hold on;
        if ~FIGUREVISIBLE, set(gcf,'Visible', 'off'); end
        hold on;
        scatter3(data(:,1),data(:,2),ones(1,numel(data(:,2)))*iTausCalculated(delayIdx),3,myColorMap(delayIdx,:),'.');%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end

    % plot
    %figure(h3); set(gcf,'Visible', FIGUREVISIBLE); clf; hold on;
    set(0, 'CurrentFigure', h3); hold on;
    offset=100; width1=500; height1=500;
    set(h3, 'Position', [(offset+width1) offset width1 height1]);     

    % scatter
    for pointIdx = 1:numel(data(:,1))
        plot(data(pointIdx,1),data(pointIdx,2),'.','Color',originColorPerTau{delayIdx}(pointIdx,:));%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end
    
    % contour (from kde)
    [bandwidth,density,X,Y] = kde2d(data);      
    [C, lCCdjk] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);
    % For later storage
    bandwidths{delayIdx}    = bandwidth;
    densities{delayIdx}     = density;
    Xs{delayIdx}            = X;
    Ys{delayIdx}            = Y;
              
    % mean
    lineH = plot(mean(data(:,1)),mean(data(:,2)),'o','MarkerFaceColor','k','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);

    title(['D# = ' num2str(iTausCalculated(delayIdx)) ', R = ' sprintf('%0.3f',correlationsPerTau(delayIdx)), 10 ,myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none')    

    xlabel(['Delta ' associatedFieldNames{1,2}] , 'Interpreter', 'none');
    ylabel(['Delta ' associatedFieldNames{1,3}] , 'Interpreter', 'none');
    %Set all fontsizes
    MW_makeplotlookbetter(15);

    xlim([  min(dataPairsPerTau{indexMidpoint}(:,1)), max(dataPairsPerTau{indexMidpoint}(:,1))  ])
    ylim([  min(dataPairsPerTau{indexMidpoint}(:,2)), max(dataPairsPerTau{indexMidpoint}(:,2))  ])

    saveas(h3,[myOutputFolder 'TIF_z_graphTauIdx_' associatedFieldNames{1,2} '_' sprintf('%05d',delayIdx) '.tif']);
    saveas(h3,[myOutputFolder 'EPS_z_graphTauIdx_' associatedFieldNames{1,2} '_' sprintf('%05d',delayIdx) '.eps'],'epsc');   

    % Let user know progress
    count=count+1;
    disp(['Finished ' num2str(count) '/' num2str(numelRangeiTausCalculated) '.'])
    
end

% Wrap the data that's needed to make this plot in a neat package :)
% ===
% Initialize
contourPlotData = struct;
% Indices of pairs calculated
contourPlotData.rangeiTausCalculated = rangeiTausCalculated;
% All available pairs, for all delays: dataPairsPerTau{tauIndex}(fieldIndex,:)
contourPlotData.dataPairsPerTau     =  dataPairsPerTau;
% Time between each frame, according to corrData (not necessarily the same, since might be extrapolation).
contourPlotData.deltaXCorrData      = deltaXCorrData;
% output of kde2d(data); plot contour lines using:
% >>contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);
contourPlotData.bandwidths           = bandwidths;
contourPlotData.densities             = densities;
contourPlotData.Xs                   = Xs;
contourPlotData.Ys                   = Ys;

% average point (used for legend too)
%{
legendLines = []; previous = 0;
for i = 1:numberOfDataFiles        
    lineH = plot(mean(data(:)),mean(data(2,:)),'o','MarkerFaceColor',distinguishableColors(myGrouping(i)+1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',7);
if myGrouping(i) ~= previous, legendLines = [legendLines lineH]; end; previous = myGrouping(i);
end    
legend( legendLines, legendDescriptions,'Location','northeast');
%}

if PLOT3DSCATTER && ~CALCULATEONLY
    %figure(h2); set(gcf,'Visible', FIGUREVISIBLE);
    set(0, 'CurrentFigure', h2); hold on;
    xlabel('Growth rate (dbl/hr)');
    ylabel('Concentration (a.u.)');
    %Set all fontsizes
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
    set(gca,'FontSize',15);
    %ylim([-750, 2000])
    %xlim([0, max([myData(:).selected_growth_rates])])
end
    
end


%% Determine raw (non-weighed) correlation per branch..
MAXLAGS=40; % width of the correlation function, i.e. frames lag 

% Plot all branches
h4=figure('Visible', FIGUREVISIBLE); clf; hold on;
numelBranches = numel(branchData);
lengthCorr = MAXLAGS*2+1;
meanR = zeros(1,lengthCorr); meanTau = zeros(1,lengthCorr); 
for branchIdx = 1:numelBranches
    
    % Gather data
    field1Data = branchData(branchIdx).(associatedFieldNames{2});
    field2Data = branchData(branchIdx).(associatedFieldNames{3});
    
    % Substract mean
    field1Data = field1Data-mean(field1Data);
    field2Data = field2Data-mean(field2Data);
    
    % Get correlations per branch
    [RThisBranch, tauThisBranch] = xcorr(...
        field2Data, ...
        field1Data, ...        
        MAXLAGS,'coeff');
        % for some reason the X and Y are inverted in xcorr function, see
        % example in documentation.
    
    l = plot(tauThisBranch, RThisBranch,'-','Color',distinguishableColors(branchIdx,:));
    %set(l, 'LineWidth', (numelBranches-branchIdx+1)/numelBranches*10);
    
    meanR   = meanR     + RThisBranch./numelBranches;
    
end

% Overlay axes
plot([-MAXLAGS,MAXLAGS] ,[0,0]  , '-k', 'LineWidth', 2);
plot([0,0]              ,[-1,1] , '-k', 'LineWidth', 2);

plot([-MAXLAGS:MAXLAGS], meanR,'k-','LineWidth',4)

xlabel('Lag (in frames)')
ylabel('Correlation')
title(['xcorr(',associatedFieldNames{2},',..',10,associatedFieldNames{3},')'],'Interpreter','None')
MW_makeplotlookbetter(20);

%% Calculate noise when an autocorrelation function is made
% check if fields 2 and 3 are identical, i.e. autocorr has been made
if strcmp(associatedFieldNames{2},associatedFieldNames{3})
    
    %% then calculate noise for this field
    rawData = [branchData(:).(associatedFieldNames{2})];
    
    % Calculate noise (std / mean)
        % noise the prefix noise is used elsewhere, but this refers to
        % mean-subtracted data.
    theStd = std(rawData);
    theMean = mean(rawData);
    theNoise = theStd/theMean;
    
    disp(['Noise calculated for ' associatedFieldNames{2}]);
    
end

%% Some summary parameters

output.numelSchnitzcells = numel(s_rm);

%% Some random checks

if PERFORMSOMECHECKS
    % distribution of times 
    % equivalent of 'hist' calculated manually because times have interval,
    % which complicates use of hist function.
    fieldcount=0; R = struct;
    % loop over timefields of interest
    for myfield = {'Y_time','MWDJK_dY_time', 'time'}

        % administration
        myfield = char(myfield);
        fieldcount = fieldcount+1;

        % get all possible values 
        R(fieldcount).YValues = unique([s_rm.(myfield)]);
        % count how many of each values are in the schnitz struct
        R(fieldcount).countYValues = [];
        for value=R(fieldcount).YValues 
            R(fieldcount).countYValues(end+1) = sum([s_rm.(myfield)]==value);
        end

    end
    % plot first two fields of interest
    figure, clf; hold on;
    if numel(R)>2, plot(R(3).YValues, R(3).countYValues, 'xk'); end
    plot(R(1).YValues, R(1).countYValues, 'or','LineWidth',2);
    plot(R(2).YValues, R(2).countYValues, 'ob','LineWidth',2);
end

%% Also create a plot of the noise over time (i.e. distribution width)
% Do this only if we're calculating an autocorrelation function (to save
% computation time)

% Note: branchData is already ignoring schnitzes to be removed and time
% outside the time period of interest, but we can't use that here (easily)
% since it contains redundant data.
%     So, we need to use s_rm and also select for the correct fitTime

if strcmp(associatedFieldNames{2},associatedFieldNames{3})
    
    fieldIndex=2;

    %%
    timeField       = associatedFieldNames{1};
    fieldOfInterest = associatedFieldNames{fieldIndex};

    % Group all data together
    allTimes     = [s_rm.(timeField)];
    allData      = [s_rm.(fieldOfInterest)];
    % make selection
    indicesToSelect = allTimes>fitTime(1) & allTimes<fitTime(2);
    theTimes     = allTimes(indicesToSelect);
    theData      = allData(indicesToSelect);
    
    % Unique values
    uniqueTimes  = sort(unique(theTimes));
    
    %
    % calculate variances per frame and store them
    varianceValuesOverTime = NaN(1,numel(uniqueTimes)); 
    stdValuesOverTime = NaN(1,numel(uniqueTimes)); 
    meanValuesOverTime = NaN(1,numel(uniqueTimes));
    numberOfValues = NaN(1,numel(uniqueTimes));
    for timeIdx = 1:numel(uniqueTimes)

        timePoint = uniqueTimes(timeIdx);

        selectionIdxs = theTimes==timePoint;    
        valuesThisFrame = theData(selectionIdxs);
        valuesThisFrame = valuesThisFrame(~isnan(valuesThisFrame));
        
        varianceValuesOverTime(timeIdx) = var(valuesThisFrame);
        stdValuesOverTime(timeIdx) = std(valuesThisFrame);
        meanValuesOverTime(timeIdx) = mean(valuesThisFrame);
        
        numberOfValues(timeIdx)  = numel(valuesThisFrame);
        
        % valuesThisFrameNormalized = valuesThisFrame./mean(valuesThisFrame);
        % variancesNormalized(timeIdx) = var(valuesThisFrameNormalized);
    end

    % calculate coefficient of variation (cv)
    coefficientOfVariationOverTime = stdValuesOverTime./meanValuesOverTime;       
    
    % calculate an estimated variance based on last 10 frames (or if total
    % <10, take all frames)
    framesToSubtract = min(numel(varianceValuesOverTime)-1,9);
    meanVarianceLast10  = mean(varianceValuesOverTime(end-framesToSubtract:end));
    medianVarianceLast10= median(varianceValuesOverTime(end-framesToSubtract:end));            
    stdVarianceLast10   = std(varianceValuesOverTime(end-framesToSubtract:end));

    % calculate an estimated cv based on last 10 frames
    meanCoefficientOfVariationLast10   = mean(coefficientOfVariationOverTime(end-framesToSubtract:end));
    medianCoefficientOfVariationLast10 = median(coefficientOfVariationOverTime(end-framesToSubtract:end));            
    stdCoefficientOfVariationLast10    = std(coefficientOfVariationOverTime(end-framesToSubtract:end));
    
    %% Create a figure
    h6=figure('Visible', FIGUREVISIBLE); clf; hold on;
    %bar(uniqueTimes/60,variancesNormalized);%,'-o','LineWidth',2)
    hCVoverTime = plot(uniqueTimes/60,coefficientOfVariationOverTime,'-','LineWidth',2);          
    title(['Population based variance' 10 'Median last 10 frms = ' num2str(medianCoefficientOfVariationLast10) ' +/- ' num2str(stdCoefficientOfVariationLast10) '.']);
    xlabel('Time (hrs)'); 
    ylabel('Coefficient of variation');
    MW_makeplotlookbetter(10);

    MYLOWN = 25;
    CVsFromSomeDataPoints = coefficientOfVariationOverTime(numberOfValues>MYLOWN);
    myYlimCV=[min([CVsFromSomeDataPoints,0]),max(CVsFromSomeDataPoints)*1.05];
    if any(numberOfValues>MYLOWN) % if there was selection use it for ylim
        ylim(myYlimCV);
    end
    
    %% Save it
    saveas(h6,[myOutputFolder 'FIG_CVovertime_' associatedFieldNames{1,fieldIndex} '.fig']);
    saveas(h6,[myOutputFolder 'TIF_CVovertime_' associatedFieldNames{1,fieldIndex} '.tif']);
    saveas(h6,[myOutputFolder 'SVG_CVovertime_' associatedFieldNames{1,fieldIndex} '.svg']);

    % Now small version
    set(hCVoverTime,'LineWidth',1);    
    saveas(h6,[myOutputFolder 'FIG_CVovertime_' associatedFieldNames{1,fieldIndex} '_small.fig']);
    saveas(h6,[myOutputFolder 'TIF_CVovertime_' associatedFieldNames{1,fieldIndex} '_small.tif']);
    saveas(h6,[myOutputFolder 'SVG_CVovertime_' associatedFieldNames{1,fieldIndex} '_small.svg']);

    %% Save data for later usage
    
    % Base info
    output.CV.(fieldOfInterest).time                   = uniqueTimes;
    output.CV.(fieldOfInterest).varianceValuesOverTime = varianceValuesOverTime;
    output.CV.(fieldOfInterest).stdValuesOverTime      = stdValuesOverTime;
    output.CV.(fieldOfInterest).meanValuesOverTime     = meanValuesOverTime;
    output.CV.(fieldOfInterest).numberOfValues         = numberOfValues;
                                                       

    % CV over time
    output.CV.(fieldOfInterest).coefficientOfVariationOverTime = coefficientOfVariationOverTime;
    
    % Variance based on last 10 frames
    output.CV.(fieldOfInterest).meanVarianceLast10   = meanVarianceLast10;
    output.CV.(fieldOfInterest).medianVarianceLast10 = medianVarianceLast10;
    output.CV.(fieldOfInterest).stdVarianceLast10    = stdVarianceLast10;

    % CV based on last 10 frames
    output.CV.(fieldOfInterest).meanCoefficientOfVariationLast10   = meanCoefficientOfVariationLast10;
    output.CV.(fieldOfInterest).medianCoefficientOfVariationLast10 = medianCoefficientOfVariationLast10
    output.CV.(fieldOfInterest).stdCoefficientOfVariationLast10    = stdCoefficientOfVariationLast10;
    
end 
   





