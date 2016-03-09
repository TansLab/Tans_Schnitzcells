% DJK_divide_branch_data divides branch data in seperate groups, each
% originating from the same parent cell. 
% Modification wrt to DJK_tim_branch_data: branches are trimmed to a
% certain initiation time "inittime".
% All branches which have the same parent cell from "DJK_trim_branchData"
% will appear in the same branchgroup. Even if inittime restricts these
% branches to a shorter time window where there is no parentschnitz in
% common.
%
% returns a collection of branches (a struct of structs)
%
function branch_groups = DJK_divide_branch_data(branches,timefield,inittime);

% FIRST CHECK HOW MANY BRANCHES IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(branches)
% comment NW 2013-04-11: if CrossCorrs will be extented to branches with
% SameLength=0, then the following lines needs to be changed (it's
% assumed that all branches have same length).
  schnitzNrs(i,:) = branches(i).schnitzNrs;
end
parent_schnitzes = schnitzNrs(:,1)';
unique_parent_schnitzes = unique(parent_schnitzes);
disp(['There are ' num2str(length(unique_parent_schnitzes)) ' unique branches in data (parent cells: ' num2str(unique_parent_schnitzes) ')']);


% investigate first branch and find idx of the time points larger than
% inittime
timevec=branches(1).(timefield);
idxtime=find(timevec>inittime);
disp(['Will restrict to initiation time ' num2str(inittime) '.']);



% DIVIDE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branch_groups = struct;
for group = 1:length(unique_parent_schnitzes)
  parent_cell = unique_parent_schnitzes(group);
  idx = find(parent_schnitzes == parent_cell);

  
  % TRIM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentbranches=branches(idx);
    trimmed_branches=currentbranches;
    fields = fieldnames(branches);
    for f = 1:length(fields)
        for i=1:length(trimmed_branches)
            trimmed_branches(i).(char(fields(f))) = currentbranches(i).(char(fields(f)))(idxtime);
        end
    end

  
  branch_groups(group).branches = trimmed_branches;
  branch_groups(group).parent_cell = parent_cell;
  branch_groups(group).nr_branches = length(idx);
end
