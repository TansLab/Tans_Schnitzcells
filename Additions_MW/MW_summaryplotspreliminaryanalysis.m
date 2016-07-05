


%% Creating overview of the data generated during preliminary analysis

if ~exist('FLUORIDXTOPLOT','var')
    error('Please set FLUORIDXTOPLOT');
    %FLUORIDXTOPLOT = [1]; % use numbers
    %FLUORIDXTOPLOT = [1 2]; % use numbers
end

% PARAMETERS ==============================================================
% Load the .mat file from the appropriate analysis.
% (summaryParametersPreliminary.mat should contain "thedata" struct with all
% "p" and "settings" parameters created during your analyses.)
if ~exist('MYDIR','var')
    error('MYDIR NOT DEFINED');
    % use 
    % MYDIR = 'F:\Datasets\2016-03-04\outputSummary\';
    % MYDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-06-21_cAMP60uM\outputSummary\';
    % to set MYDIR.
end
load([MYDIR 'summaryParametersPreliminary.mat'])

% optional input
% MANUALSCHNITZLOCATIONS = {'..','..'};

% END PARAMETERS ==========================================================

% Load data per position
clear alldata;
numelThedata = numel(thedata); positionIndices = [];
for i = 1:numelThedata

    if isempty(thedata(i).summaryParameters)
        continue
    else
        positionIndices(end+1) = i;
    end
    
    currentP = thedata(i).p;
    % ugly work around for when stored schnitz location is not correct any more
    if exist('MANUALSCHNITZLOCATIONS','var')
        currentP.schnitzName = MANUALSCHNITZLOCATIONS{i};
        warning('Ugly solution for schnitz location issuel;');
    end
    [currentP,currentSchnitzcells] = DJK_compileSchnitzImproved_3colors(currentP,'quickMode',1);
    
    % create struct with summary data per frame
    % - frameData(i).frame_nr gives frame_nr for record i.
    % - frameData(i).<data> gives data for record i with framenumber frameData(i).frame_nr.
    frameData = struct; count=0;
    for frameNr = unique([currentSchnitzcells.frame_nrs])
        
        % Administration
        count=count+1;        
        frameData(count).frame_nr = frameNr;
        
        % Select relevant data for this frame
        %{
        %applicable to preliminary proxy schnitzcells structs only
        indexesForThisFrame = find([currentSchnitzcells(:).frame_nrs] == frameNr);
        currentFrameSchnitzes = currentSchnitzcells(indexesForThisFrame);
        %}
        % applicable both preliminary and full analysis schnitzcells
        % structs
        % Loop over schnitzes and see whether it contains records that
        % belong to the current frame of interest.
        disp(['Calculating frame ' num2str(frameNr)]);
        schnitzNrAndIndexForCurrentFrame=[]; %[schnitznr,idx;..]
        relevantSchnitzes = []; hitCount=0;
        allFields = fieldnames(currentSchnitzcells(1)); numelallFields = numel(allFields);
        for schnitzIdx=1:numel(currentSchnitzcells)
            
            frameMatch = find(currentSchnitzcells(schnitzIdx).frame_nrs==frameNr);
            
            % record relevant indices when there's a match
            if ~isempty(frameMatch)
                hitCount=hitCount+1;
                schnitzNrAndIndexForCurrentFrame=[schnitzNrAndIndexForCurrentFrame;schnitzIdx,frameMatch];
            
                % record relevant schnitzinfo at that framenumber            
                for schnitzFieldIdx = 1:numelallFields
                    currentFieldName = allFields{schnitzFieldIdx};
                    relevantSchnitzes(hitCount).(currentFieldName) = ...
                            currentSchnitzcells(schnitzIdx).(currentFieldName);
                end
            end
        end
        
        % Create summary parameters for this frame
        % ===
        
        % Time at which frame was recorded
        schnitzLifeTimeIdx=schnitzNrAndIndexForCurrentFrame(1,2); % which record to take for first frame
        frameData(count).frameTime = relevantSchnitzes(1).time(schnitzLifeTimeIdx); % all times for frame should be equal, select 1st
        
        % Summed colony length
        frameSummedLength=0;
        for idx=1:hitCount
            frameSummedLength = frameSummedLength + ...
                relevantSchnitzes(schnitzNrAndIndexForCurrentFrame(idx,1)).length_fitNew(schnitzNrAndIndexForCurrentFrame(idx,2));
        end
        frameData(count).frameSummedLength = frameSummedLength
        %frameData(count).frameSummedLength = sum([relevantSchnitzes.length_fitNew]);      
        
        % Colony fluor mean
        for fluorIdx = FLUORIDXTOPLOT            
            fluorFieldName = [upper(currentP.(['fluor' num2str(fluorIdx)])(1)) '5_mean_all'];
            
            frameFluorSum=0;
            for idx=1:hitCount
                frameFluorSum = frameFluorSum + ...
                    relevantSchnitzes(schnitzNrAndIndexForCurrentFrame(idx,1)).(fluorFieldName)(schnitzNrAndIndexForCurrentFrame(idx,2));
            end
            frameData(count).frameFluorSum(fluorIdx) = frameFluorSum;
            %frameData(count).fluorMean(fluorIdx) = mean([relevantSchnitzes.(fluorFieldName)]);
        end
                
    end
    
    % Save schnitzcells info per frame for position i
    alldata(i).frameData = frameData;
    alldata(i).label = [currentP.movieDate currentP.movieName];
    
    disp(['Finished schnitzcells ' num2str(i) '/' num2str(numelThedata) '']);
    
end

%% Prepare plotting
FONTSIZE=15;
% MAKE CUSTOMCOLORS BELOW MATCH POSITIONS!

mycolors=colormap('lines');
preferredcolors=[0 0 0;0 0.568627450980392 1;0 0.709803921568627 0.235294117647059;0.545098039215686 0.450980392156863 0.92156862745098;0.564705882352941 0.552941176470588 0.619607843137255;0.909803921568627 0 0];% from some_colors.m
preferredcolors=[preferredcolors;preferredcolors;preferredcolors];
if ~exist('CUSTOMCOLORS','var')
    CUSTOMCOLORS= [preferredcolors(2,:);preferredcolors(2,:);preferredcolors(2,:);...
                   preferredcolors(3,:);preferredcolors(3,:);preferredcolors(3,:)]; % make this match positions
end


%% Go over positions again and now create fluor plots======================

figure(1); clf;

plotrows=numel(FLUORIDXTOPLOT);

fluorValuesAll=[];
for i = positionIndices
for frIdx = 1:numel(alldata(i).frameData)
for fluorIdx = FLUORIDXTOPLOT
    fluorValuesAll{i}(fluorIdx,frIdx) = alldata(i).frameData(frIdx).fluorMean(fluorIdx);
end
end
end

fluorYlim1={}; 
for fluorIdx = FLUORIDXTOPLOT

    subplot(plotrows,2,1+(2*(fluorIdx-1)) ); hold on;

    % Plot raw fluor data------------------------------------------------------
    for i = positionIndices
        fluorValuesCurrent = squeeze(fluorValuesAll{i}(fluorIdx,:));
        l = plot([alldata(i).frameData(:).frameTime],fluorValuesCurrent,'-x');
        set(l, 'LineWidth', 2, 'Color', CUSTOMCOLORS(i,:));

    end

    % Determine ylim & xlim
    fluordatapile=[];timeDataPile=[];
    for i=positionIndices
        fluorValuesCurrent = squeeze(fluorValuesAll{i}(fluorIdx,:));
        fluordatapile = [fluordatapile, fluorValuesCurrent];
        timeDataPile = [timeDataPile, alldata(i).frameData(:).frameTime];
    end
    maxtimeDataPile=max(timeDataPile);
    fluorYlim1{fluorIdx} = [0,max(fluordatapile)*1.1];
    ylim(fluorYlim1{fluorIdx});
    xlim([0,maxtimeDataPile]);

    % Set title etc.
    title('')
    MW_makeplotlookbetter(FONTSIZE);
    xlabel('time (min)');
    ylabel('signal (a.u.)');

    % Plot normalized fluor data-----------------------------------------------
    subplot(plotrows,2,2+(2*(fluorIdx-1)) ); hold on;

    lines = [];
    for i = positionIndices

        % get median
        fluorValuesCurrent = squeeze(fluorValuesAll{i}(fluorIdx,:))';
        %fluorValuesCurrent = [alldata(i).frameData(:).fluorMean];
        notNanIndices = ~isnan(fluorValuesCurrent);
        fluorCurrentmedian = median(fluorValuesCurrent(notNanIndices));

        % current Y data to plot
        currentFluorY = fluorValuesCurrent./fluorCurrentmedian;    
        sizecurrentFluorY = size(currentFluorY);
        
        % transpose data if required, just for plotting
        % (this is a bit ugly, but required for backward compatibility)
        if sizecurrentFluorY(1)>sizecurrentFluorY(2)
            currentFluorY=currentFluorY';
        end
            
        % plot data
        currentTimes = [alldata(i).frameData(:).frameTime];
        l = plot(currentTimes(notNanIndices),currentFluorY(notNanIndices),'-x');
        % cosmetics
        set(l, 'LineWidth', 2, 'Color', mycolors(i,:));
        lines(end+1)=l; % for legend later

    end
    xlim([0,maxtimeDataPile]);
    fluorYlim2 = [0,2];
    ylim(fluorYlim2);

    % Line at 1 (lines are median-normalized, should be ~1)
    plot([0,maxtimeDataPile],[1,1],'-','Color','k','LineWidth',3)

    % Set title etc.
    title('')
    MW_makeplotlookbetter(FONTSIZE);
    xlabel('time (min)');
    ylabel('signal (median-normalized)');
    legend(lines, {alldata(positionIndices).label},'Location', 'NorthOutside');

    % save the figure 
    saveas(1, [MYDIR 'summaryPlot_fluorsignal.tif'],'tif');
    saveas(1, [MYDIR 'summaryPlot_fluorsignal.eps'],'epsc');
end

%% Go over positions again and now create length/growth plots==============
figure(2); clf; 

flag=0;
if exist('settings','var') 
    if isfield(settings,'fitTimeMu')
        FITWINDOW = settings.fitTimeMu;
        flag=1;
    end
end
if ~flag
    FITWINDOW = [250,550];
    warning('FITWINDOW set to default');
end

subplot(1,2,1); hold on;

% Plot summed length data--------------------------------------------------
fittedMus = zeros(1,max(positionIndices));
for i = positionIndices
        
    % the data for all times
    time = [alldata(i).frameData(:).frameTime]./60; % convert to hours
    datatofit = log([alldata(i).frameData(:).frameSummedLength])./log(2); % fit exponential w. base 2
    % determine which data to fit (i.e. within fitting window)
    fitWindowHrs = FITWINDOW./60;
    indexesInTimeWindow = time>fitWindowHrs(1) & time<fitWindowHrs(2);
    % fit it
    myfit = polyfit(time(indexesInTimeWindow),datatofit(indexesInTimeWindow),1);
    fittedMus(i)=myfit(1)
    
    % plot it
    l = plot([alldata(i).frameData(:).frameTime],[alldata(i).frameData(:).frameSummedLength],'-x');
    set(l, 'LineWidth', 2, 'Color', CUSTOMCOLORS(i,:));        
end

set(gca,'yscale','log');

% Determine ylim
frameSummedLengthdatapile=[];
for i=positionIndices
    frameSummedLengthdatapile = [frameSummedLengthdatapile, alldata(i).frameData(:).frameSummedLength];
end
theYLim = [min(frameSummedLengthdatapile)/2,max(frameSummedLengthdatapile)*2];

% Plot timewindow used for fitting
plot(FITWINDOW,[theYLim(1) theYLim(1)],'k-^','LineWidth',2,'MarkerFaceColor','k');
%plot(FITWINDOW(1),theYLim(1),'k>','LineWidth',2,'MarkerFaceColor','k')
%plot(FITWINDOW(2),theYLim(1),'k<','LineWidth',2,'MarkerFaceColor','k')

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('summed length (a.u.)');
ylim(theYLim);

% bar plot of fitted growth rates------------------------------------------
subplot(1,2,2); hold on;

% plot
for i=positionIndices
    l=bar(i, fittedMus(i));
    set(l, 'LineWidth', 1, 'FaceColor', CUSTOMCOLORS(i,:), 'EdgeColor', CUSTOMCOLORS(i,:));
end

% align,labeling
set(gca, 'XTickLabel',{1:max(positionIndices)}, 'XTick',1:max(positionIndices))
xlim([0,max(positionIndices)+1])

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('fitted growth rate (doublings/hour)');

% save the figure 
saveas(2, [MYDIR 'summaryPlot_colonylength.tif'],'tif');
saveas(2, [MYDIR 'summaryPlot_colonylength.eps'],'epsc');

%% plot for 2 growth regimes







