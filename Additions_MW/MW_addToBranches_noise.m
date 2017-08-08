function branchData = MW_addToBranches_noise(branchData, associatedFieldNames)
% function branchData = MW_addToBranches_noise(branchData, associatedFieldNames)
% 
% Input: 
%   branchData:             branch structure
%   associatedFieldNames:   cell with names of timefield, fieldX, fieldY
%                           respectively
%
% Output:
%   relative noise field, which is original field with prefix "noise_".
%
% Note that this output is different from the DJK_addToBranches_noise
% function, which outputs three different types of normalization. The
% advantage of this function though is that it can handle branches of
% unequal width.


allTimePoints = unique([branchData(:).(associatedFieldNames{1})]);

    for timeIdx = 1:numel(allTimePoints) 

        %%

        % gather the data for this timepoint

        % Current timepoint
        currentTimePoint=allTimePoints(timeIdx);    

        % Collect indices for each branch that match with this timepoints
        % (empty if timepoint does not exist in this data)
        timeHits= arrayfun(@(br) find(branchData(br).(associatedFieldNames{1}) == currentTimePoint), 1:numel(branchData),'UniformOutput',0);
        % Branches for which there is data
        hitBranches= find(arrayfun(@(br) ~isempty(find(branchData(br).(associatedFieldNames{1}) == currentTimePoint)), 1:numel(branchData)));

        % get the data
        data2ForThisTimePoint = arrayfun(@(br) branchData(br).(associatedFieldNames{2})(timeHits{br}), hitBranches);
        data3ForThisTimePoint = arrayfun(@(br) branchData(br).(associatedFieldNames{3})(timeHits{br}), hitBranches);

        % get the mean
        data2Mean = mean(data2ForThisTimePoint);
        data3Mean = mean(data3ForThisTimePoint);

        % put the data back
        for brIdx = 1:numel(hitBranches)
            br=hitBranches(brIdx);
            branchData(br).(['noise_' associatedFieldNames{2}])(timeHits{br}) = data2ForThisTimePoint(brIdx)./data2Mean-1;
            branchData(br).(['noise_' associatedFieldNames{3}])(timeHits{br}) = data3ForThisTimePoint(brIdx)./data3Mean-1;
        end

        %{
        % get applicable indices for this point in time
        currentTimePoint=allTimePoints(currentTimePointIdx);


        % Now we know which branches to mix, so do that
        randomizedBranches = hitBranches(randperm(numel(hitBranches)));

        for brIdx=1:numel(hitBranches)

            sourceBranchIdx = hitBranches(brIdx);
            targetBranchIdx = randomizedBranches(brIdx);

            % actual randomization            
            branch_groupsControl(1).branches(sourceBranchIdx).([associatedFieldNames{3} CONTROLSUFFIX])(timeHits{sourceBranchIdx}) = ...
                branch_groups(1).branches(targetBranchIdx).([associatedFieldNames{3}])(timeHits{targetBranchIdx});
            % also do for normalized field
            branch_groupsControl(1).branches(sourceBranchIdx).([FIELDPREFIX associatedFieldNames{3} CONTROLSUFFIX])(timeHits{sourceBranchIdx}) = ...
                branch_groups(1).branches(targetBranchIdx).([FIELDPREFIX associatedFieldNames{3}])(timeHits{targetBranchIdx});

        end
        %}
    end
    
end