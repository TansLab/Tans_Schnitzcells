


%% Creating overview of the data generated during preliminary analysis

% PARAMETERS ==============================================================
% Load the .mat file from the appropriate analysis.
% (summaryParametersPreliminary.mat should contain "thedata" struct with all
% "p" and "settings" parameters created during your analyses.)
MYDIR = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-12-02\outputSummary\';
load([MYDIR 'summaryParametersPreliminary.mat'])

fluorsIdxToPlot = [1]; % use capital letter abbreviations
% END PARAMETERS ==========================================================

% Load data per position
numelThedata = numel(thedata);
for i = 1:numelThedata
    
    [currentP,currentSchnitzcells] = DJK_compileSchnitzImproved_3colors(thedata(i).p,'quickMode',1);
    
    % create struct with summary data per frame
    % - frameData(i).frame_nr gives frame_nr for record i.
    % - frameData(i).<data> gives data for record i with framenumber
    %     frameData(i).frame_nr.
    frameData = struct; count=0;
    for frameNr = unique([currentSchnitzcells.frame_nrs])
        
        % Administration
        count=count+1;        
        frameData(count).frame_nr = frameNr;
        
        % Select relevant data for this frame
        indexesForThisFrame = find([currentSchnitzcells(:).frame_nrs] == frameNr);
        currentFrameSchnitzes = currentSchnitzcells(indexesForThisFrame);
        
        % Create summary parameters for this frame
        % ===
        
        % Time at which frame was recorded
        frameData(count).frameTime = currentFrameSchnitzes(1).time; % all times for frame should be equal, select 1st
        % Summed colony length
        frameData(count).frameSummedLength = sum([currentFrameSchnitzes.length_fitNew]);
        
        % Colony fluor mean
        for fluorIdx = fluorsIdxToPlot            
            fluorFieldName = [upper(currentP.(['fluor' num2str(fluorIdx)])(1)) '5_mean_all'];
            frameData(count).fluorMean(fluorIdx) = mean([currentFrameSchnitzes.(fluorFieldName)]);
        end
                
    end
    
    % Save schnitzcells info per frame for position i
    alldata(i).frameData = frameData;
    alldata(i).label = [currentP.movieDate currentP.movieName];
    
    disp(['Finished schnitzcells ' num2str(i) '/' num2str(numelThedata) '']);
    
end

%% Plot
FONTSIZE=15;

mycolors=colormap('lines');
preferredcolors=[0 0 0;0 0.568627450980392 1;0 0.709803921568627 0.235294117647059;0.545098039215686 0.450980392156863 0.92156862745098;0.564705882352941 0.552941176470588 0.619607843137255;0.909803921568627 0 0];% from some_colors.m
preferredcolors=[preferredcolors;preferredcolors;preferredcolors];
customcolors= [preferredcolors(2,:);preferredcolors(2,:);preferredcolors(2,:);preferredcolors(2,:);...
               preferredcolors(3,:);preferredcolors(3,:);preferredcolors(3,:);preferredcolors(3,:)]; % make this match positions


%% Go over positions again and now create fluor plots======================
figure(1); clf; 

subplot(1,2,1); hold on;

% Plot raw fluor data------------------------------------------------------
for i = 1:numelThedata
    
    l = plot([alldata(i).frameData(:).frameTime],[alldata(i).frameData(:).fluorMean],'-');
    set(l, 'LineWidth', 4, 'Color', customcolors(i,:));

end

% Determine ylim & xlim
fluordatapile=[];timeDataPile=[];
for i=1:numel(alldata)
    fluordatapile = [fluordatapile, alldata(i).frameData(:).fluorMean];
    timeDataPile = [timeDataPile, alldata(i).frameData(:).frameTime];
end
maxtimeDataPile=max(timeDataPile);
ylim([0,max(fluordatapile)*1.1]);
xlim([0,maxtimeDataPile]);

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('signal (a.u.)');

% Plot normalized fluor data-----------------------------------------------
subplot(1,2,2); hold on;

for i = 1:numelThedata
    
    % get median
    fluorCurrentmedian = median([alldata(i).frameData(:).fluorMean]);
    
    % current Y data to plot
    currentFluorY = [alldata(i).frameData(:).fluorMean]./fluorCurrentmedian;    
    
    % plot data
    l = plot([alldata(i).frameData(:).frameTime],currentFluorY,'-');
    % cosmetics
    set(l, 'LineWidth', 3, 'Color', mycolors(i,:));
    
end
xlim([0,maxtimeDataPile]);
ylim([0,2]);

% Line at 1 (lines are median-normalized, should be ~1)
plot([0,maxtimeDataPile],[1,1],'-','Color','k','LineWidth',3)

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('signal (median-normalized)');


% save the figure 
saveas(1, [MYDIR 'summaryPlot_fluorsignal.tif'],'tif');
saveas(1, [MYDIR 'summaryPlot_fluorsignal.eps'],'epsc');

%% Go over positions again and now create fluor plots======================
figure(2); clf; 

subplot(1,2,1); hold on;

% Plot raw fluor data------------------------------------------------------
for i = 1:numelThedata    
    
    % fit mu
    time = [alldata(i).frameData(:).frameTime]./60; % convert to hours
    datatofit = log([alldata(i).frameData(:).frameSummedLength])./log(2); % fit exponential w. base 2
    myfit = polyfit(time,datatofit,1);
    fittedMus(i)=myfit(1)
    
    l = plot([alldata(i).frameData(:).frameTime],[alldata(i).frameData(:).frameSummedLength],'-');
    set(l, 'LineWidth', 4, 'Color', customcolors(i,:));

end

set(gca,'yscale','log');

% Determine ylim
frameSummedLengthdatapile=[];
for i=1:numel(alldata)
    frameSummedLengthdatapile = [frameSummedLengthdatapile, alldata(i).frameData(:).frameSummedLength];
end
ylim([0,max(frameSummedLengthdatapile)*1.1]);

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('summed length (a.u.)');

% Plot normalized fluor data-----------------------------------------------
subplot(1,2,2); hold on;

% plot
for i=1:numel(alldata)
    l=bar(i, fittedMus(i));
    set(l, 'LineWidth', 1, 'FaceColor', customcolors(i,:), 'EdgeColor', customcolors(i,:));
end

% align,labeling
set(gca, 'XTickLabel',{1:numel(alldata)}, 'XTick',1:numel(alldata))
xlim([0,numel(alldata)+1])

% Set title etc.
title('')
MW_makeplotlookbetter(FONTSIZE);
xlabel('time (min)');
ylabel('fitted growth rate (doublings/hour)');

% save the figure 
saveas(2, [MYDIR 'summaryPlot_colonylength.tif'],'tif');
saveas(2, [MYDIR 'summaryPlot_colonylength.eps'],'epsc');









