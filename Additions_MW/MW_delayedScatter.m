
%% Description

% MW 2015/04
%
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
%
% Additional (optional) parameters:
% - NOTMEANSUBTRACTED   : Don't automatically take fields that have the
%                         noise subtracted, i.e. noise_FieldOfInterest.
% - FIGUREVISIBLE       : hides some figures, only saves them 
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

PERFORMSOMECHECKS = 0;

if ~exist('myOutputFolder')
    myOutputFolder = 'C:\Users\wehrens\Desktop\testdelayedscatter\output\';
end

if ~exist('associatedFieldNames') | ~exist('p') | ~exist('badSchnitzes')
    error('input not supplied.')
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
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_Conc_oldRates'];


%% Plot branches
HIGHLIGHTSUSPICOUS = 0;
%yfieldbranchtoplot = 3; % 3 = growth, 2 = fluor

for yfieldbranchtoplot=[2,3]

    % Just some plot colors
    distinguishableColors = distinguishable_colors(numel(branchData)+1,[1 1 1]); 

    % Plot all branches
    h1=figure(); clf; hold on;
    offset=100; width1=800; height1=600;
    set(h1, 'Position', [offset offset width1 height1]);
    numelBranches = numel(branchData);
    for branchIdx = 1:numelBranches
        l = plot(branchData(branchIdx).(associatedFieldNames{1}), branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot}),'-o','Color',distinguishableColors(branchIdx,:))
        set(l, 'LineWidth', (numelBranches-branchIdx+1)/numelBranches*10);
    end

    % Brute force average branches (note weighing is irrelevant for average)
    branchMatrixFieldX = []; branchMatrixFieldY = [];
    for branchIdx = 1:numelBranches
        currentXvector = branchData(branchIdx).(associatedFieldNames{1});
        branchMatrixFieldX = [branchMatrixFieldX; currentXvector];

        currentYvector = branchData(branchIdx).(associatedFieldNames{yfieldbranchtoplot});
        branchMatrixFieldY = [branchMatrixFieldY; currentYvector];
    end
    meanXvector = mean(branchMatrixFieldX);
    meanYvector = mean(branchMatrixFieldY);

    % xlabel
    xlabel(associatedFieldNames{1},'Interpreter', 'None'), ylabel(associatedFieldNames{yfieldbranchtoplot},'Interpreter', 'None')

    myXlimFig1 = max(branchData(branchIdx).(associatedFieldNames{1}));
    xlim([0, myXlimFig1]);
    myYlimFig1 = [min([0, [branchData.(associatedFieldNames{yfieldbranchtoplot})]]),...
                  max([branchData.(associatedFieldNames{yfieldbranchtoplot})])];
    ylim([myYlimFig1(1), myYlimFig1(2)*1.5]);


    %Set all fontsizes
    MW_makeplotlookbetter(20);
    
    % Plot histogram
    h2=figure(), clf, hold on
    allYdata = [branchData.(associatedFieldNames{yfieldbranchtoplot})];
    [nelements, centers] = hist(allYdata,200)
    deltaY = centers(2)-centers(1);
    totalCount = numel(allYdata);
    %nelements=nelements./deltaY;
    %plot(centers,nelements,'or','LineWidth',2)
    bar(centers,nelements,'FaceColor','r','EdgeColor','r')
    % Fit distribution
    pd=fitdist(allYdata', 'Normal')
    fittedDistrX = [min(allYdata):(max(allYdata)-min(allYdata))/100:max(allYdata)]
    fittedDistrYnorm = normpdf(fittedDistrX,pd.mu,pd.sigma)
    fittedDistrY = fittedDistrYnorm.*totalCount.*deltaY;
    plot(fittedDistrX,fittedDistrY,'k', 'LineWidth', 3);
    %probplot(allYdata)
    MW_makeplotlookbetter(20);

    title(['PDF for ' associatedFieldNames{yfieldbranchtoplot}],'Interpreter','None');
    xlabel(associatedFieldNames{yfieldbranchtoplot},'Interpreter','None');
    ylabel('PDF(s) * N * \Deltas (counts)');

    % Define some axes limits
    myYlim = max(nelements)*1.1;
    ylim([0,myYlim]);
    xlim([myYlimFig1(1), myYlimFig1(2)*1.5]);

    % TODO make 99% confidence and plot in previous figure.
    %confidence = paramci(pd,'Alpha',.01);
    sigma2 = pd.mu + 4.*[-pd.sigma, pd.sigma];
    plot([sigma2(1),sigma2(1)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    plot([sigma2(2),sigma2(2)],[0,myYlim],'--','Color',[.5 .5 .5],'LineWidth', 2)
    sigma5 = pd.mu + 5.*[-pd.sigma, pd.sigma];
    plot([sigma5(1),sigma5(1)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)
    plot([sigma5(2),sigma5(2)],[0,myYlim],':','Color',[.5 .5 .5],'LineWidth', 2)

    % Now also plot confidence intervals in previous figure
    figure(h1), hold on;
    l1=plot([0,myXlimFig1],[sigma2(1),sigma2(1)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    plot([0,myXlimFig1],[sigma2(2),sigma2(2)],'--','Color',[.5 .5 .5],'LineWidth', 2)
    l2=plot([0,myXlimFig1],[sigma5(1),sigma5(1)],':','Color',[.5 .5 .5],'LineWidth', 2)
    plot([0,myXlimFig1],[sigma5(2),sigma5(2)],':','Color',[.5 .5 .5],'LineWidth', 2)

    legend([l1,l2],{'2\sigma confidence','5\sigma confidence'},'location','Best');

    % Now list schnitzes that have suspiciously high signal:
    mySigma = sigma2;
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

    % Plot mean behavior
    figure(h1), hold on;
    plot(meanXvector, meanYvector,'-','Color','k','LineWidth',3)
    
    saveas(h1,[myOutputFolder 'TIF_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.tif']);
    saveas(h1,[myOutputFolder 'EPS_branches_' associatedFieldNames{1,yfieldbranchtoplot} '.eps'],'epsc');

    saveas(h2,[myOutputFolder 'TIF_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.tif']);
    saveas(h2,[myOutputFolder 'EPS_PDF_' associatedFieldNames{1,yfieldbranchtoplot} '.eps'],'epsc');
    
    % For later output
    output.branchavg.([associatedFieldNames{yfieldbranchtoplot} '_xfield']) = meanXvector; % usually time
    output.branchavg.(associatedFieldNames{yfieldbranchtoplot}) = meanYvector;
    
end

%% Get the actual cross-corrs

%REDUNDANCYALLOWED = 2^2;
REDUNDANCYALLOWED = 2^2;
ONSCREEN=1;
NRBRANCHGROUPS=4;
FIELDPREFIX = 'noise_';
CONTROLSUFFIX = '_randomizedlineages';
%FIELDPREFIX = 'relative_';

% Some additional editing of the branches:
branchData = DJK_addToBranches_noise(p, branchData,'dataFields',{associatedFieldNames{1},associatedFieldNames{2},associatedFieldNames{3}});

% Trim of starting frames until there are N start schnitzes
trimmed_branches = DJK_trim_branch_data(branchData,NRBRANCHGROUPS);
% Divide branchdata in groups based on those N start schnitzes
branch_groups = DJK_divide_branch_data(trimmed_branches);

% Colony average mean has already been substracted so theoretically extra
% normalization shouldn't have an effect.
p.extraNorm=0;

% THIS MIGHT FAIL BECAUSE TIME IS NOT SET CORRECTLY (SHOULD BE time_at_..)
% To calculate cross correlations, additional normalization is usually
% performed, namely to filter out colony average behavior.
[CorrData,composite_corr] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, [FIELDPREFIX associatedFieldNames{1,2}],[FIELDPREFIX associatedFieldNames{1,3}] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 

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

[CorrDataControl, composite_corrControl] = DJK_plot_crosscorrelation_standard_error_store(p, branch_groupsControl, [FIELDPREFIX associatedFieldNames{1,2}],[FIELDPREFIX associatedFieldNames{1,3} CONTROLSUFFIX] ,'selectionName',name_rm_branch,'timeField',associatedFieldNames{1},'onScreen',ONSCREEN); 

% Do we want to filter out colony average behavior for the "delayed
% scatter" plots also? Maybe do this with noise fields?
% But let's try with "raw" data first..
p.timeField = associatedFieldNames{1,1};
% p.tauIndices = [-7:7]; %p.tauIndices = [-29:4:-1,0,1:4:30];
if isfield(p,'tauIndices'), p=rmfield(p,'tauIndices'); end
[dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = ...
    MW_getdelayedscatter(p, branchData, [FIELDPREFIX associatedFieldNames{1,2}], [FIELDPREFIX associatedFieldNames{1,3}], REDUNDANCYALLOWED)

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
    h98=figure(), clf, hold on;
    
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

myfig=figure(),clf,hold on;
l=plot(iTausCalculated,correlationsPerTau,'o-r','LineWidth',2)

%% Compare two cross-corrs (DJK & MW)

h5=figure(),clf,hold on;

% Plot DJK cross correlation function
errorbar(CorrData(:,1),CorrData(:,2),CorrData(:,3),'s-','Color', [.5,.5,.5], 'LineWidth',2)
l1=plot(CorrData(:,1),CorrData(:,2),'s-k','LineWidth',2)

% Plot control DJK cross correlation function
%errorbar(CorrDataControl(:,1),CorrDataControl(:,2),CorrDataControl(:,3),'x-','Color', [.5,.5,.5], 'LineWidth',2)
l3=plot(CorrDataControl(:,1),CorrDataControl(:,2),'x-b','LineWidth',2)

% Calculate appropriate x-axis assuming assuming same delta(x) as CorrData,
% and dx is same everywhere.
centerIdx=ceil(size(CorrData,1)/2);
deltaXCorrData = CorrData(centerIdx+1,1)-CorrData(centerIdx,1);
MWxAxis = iTausCalculated.*deltaXCorrData;

% Plot MW cross correlation function
l2=plot(MWxAxis,correlationsPerTau,'o-r','LineWidth',2)

% If you recalculate correlations again w. different params, this allows
% plotting of extra line.
%l3=plot(CorrData(:,1),correlationsPerTau100,'o-','Color',[1 .5 0],'LineWidth',2)

myxlimvalues=[min(CorrData(:,1)), max(CorrData(:,1))];
xlim(myxlimvalues);
ylim([-1.1,1.1]);
plot(myxlimvalues,[0,0],'k-');

%legend([l1,l2,l3],{'DJK','MW','Control'})
legend([l1,l2,l3],{'Lineages','Scatters','Control'})

title(['DJK vs. MW R -- ' myID '_' p. movieDate  '_' p.movieName], 'Interpreter', 'none');
xlabel('\tau (hrs)');
ylabel(['R(' associatedFieldNames{1,2} ',... ' 10 associatedFieldNames{1,3} ') (normalized)'], 'Interpreter', 'none');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal');
set(gca,'FontSize',15);

plot([0,0],[-1,1],'-k');

saveas(h5,[myOutputFolder 'TIF_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '.tif']);
saveas(h5,[myOutputFolder 'EPS_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '.eps'],'epsc');



%% Plot code from CRPcAMP..overview..general
% ==========
if ~exist('PLOTSCATTER','var'), PLOTSCATTER=1; end;
if PLOTSCATTER
NRCONTOURLINES = 5;
SHOWPLUSMINFROMZERO = 25;
PLOT3DSCATTER = 0;

% Whether plots should be shown.
if ~exist('FIGUREVISIBLE','var'), FIGUREVISIBLE=1; end;

% What range should be plotted.
numeliTausCalculated=numel(iTausCalculated);
indexMidpoint = ceil(numeliTausCalculated/2)
ShowPlusMinFromZero=min(SHOWPLUSMINFROMZERO, indexMidpoint-1);
rangeiTausCalculated = [indexMidpoint-ShowPlusMinFromZero:indexMidpoint+ShowPlusMinFromZero];
numelRangeiTausCalculated = numel(rangeiTausCalculated);
% delayIdx = 11; % 39 is middle

myColorMap = colormap(winter(numel(iTausCalculated)));

if PLOT3DSCATTER
    h2 = figure(); clf; hold on;
    offset=100; width1=500; height1=500;
    set(h2, 'Position', [offset offset width1 height1]);
end

h3 = figure(); 
if ~FIGUREVISIBLE, set(gcf,'Visible', 'off'); end
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
        h2 = figure(h2);
        if ~FIGUREVISIBLE, set(gcf,'Visible', 'off'); end
        hold on;
        scatter3(data(:,1),data(:,2),ones(1,numel(data(:,2)))*iTausCalculated(delayIdx),3,myColorMap(delayIdx,:),'.');%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end

    % plot
    figure(h3); clf, hold on;
    offset=100; width1=500; height1=500;
    set(h3, 'Position', [(offset+width1) offset width1 height1]);     

    % scatter
    for pointIdx = 1:numel(data(:,1))
        plot(data(pointIdx,1),data(pointIdx,2),'.','Color',originColorPerTau{delayIdx}(pointIdx,:));%,'Color',distinguishableColors(myGrouping(i)+1,:),'MarkerSize',3);
    end
    
    % contour (from kde)
    [bandwidth,density,X,Y] = kde2d(data);      
    [C, l1] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);
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
    figure(h2);
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
h4=figure(); clf; hold on;
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
        field1Data, ...
        field2Data, ...
        MAXLAGS,'coeff');
    
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
    
    % then calculate noise for this field
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








