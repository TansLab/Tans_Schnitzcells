



TFIELDS = {'time',              'time_atdY',     'time_atY', 'time_atdC',        'time_atC'};
%YFIELDS = {'muP5_fitNew_all',   'dY5_cycCor',    'Y6_mean',  'dC5_cycCor',       'C6_mean'};
YFIELDS = {'muP5_fitNew_all',   'dY5_cycCor',    'Y6_mean',  'dC5_cyccor',       'C6_mean'};
YLIMS   = {[-1,3],              [],              [],         [],                 []};
YLABELS = {'Time (hrs)',        'CRP reporter production','CRP reporter concentration','Constant production rate','Constant reporter concentration'};

for fieldIdx=1:numel(TFIELDS)
    currentTField = TFIELDS{fieldIdx};
    currentYField = YFIELDS{fieldIdx};
    currentYlim   = YLIMS{fieldIdx};

    figure(fieldIdx); clf; hold on;
    numelCountsY = [];
    dataCollectedT = [];
    dataCollectedY = [];
    for schnitzIdx = 1:numel(schnitzcells)

        tData = schnitzcells(schnitzIdx).(currentTField)/60;
        yData = schnitzcells(schnitzIdx).(currentYField);

        if numel(tData)==numel(yData)
            plot(tData, yData,'.','Color',[.5 .5 .5]);
            numelCountsY(end+1)=numel(yData);

            dataCollectedT = [dataCollectedT tData];
            dataCollectedY = [dataCollectedY yData];
        else
            disp(['Schnitz ' num2str(schnitzIdx) ' did not have equally sized t and y values.']);
        end

    end

    %%
    % Calculate average per timepoint
    [meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins] = binnedaveraging({dataCollectedT},{dataCollectedY},min(dataCollectedT):1/60:max(dataCollectedT));


    %% plot
    plot(binCenters,meanValuesForBins,'-ok','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k');

    % ylim
    if ~isempty(currentYlim)
            ylim(currentYlim); 
    else
        sortedYData = sort(dataCollectedY);
        numelY = numel(sortedYData);
        autoYlims = [sortedYData(round(0.01*numelY)), sortedYData(round(0.99*numelY))];
        ylim(autoYlims);
    end
    
    ylabel(YLABELS{fieldIdx});
    xlabel('Time (hrs)');
    MW_makeplotlookbetter(15);
    
end;

%xlim([]);
%ylim([0,2]);