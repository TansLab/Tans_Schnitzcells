function branch_groups = DJK_divide_branch_data(branches);
% DJK_divide_branch_data divides branch data in seperate groups, each
% originating from the same parent cell. 
%
% This is done to be able to finally obtain an estimate of the error. See
% also README.txt.
%
% returns a collection of branches (a struct of structs)
%

%% FIRST CHECK HOW MANY BRANCHES IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
for i=1:length(branches)
% comment NW 2013-04-11: if CrossCorrs will be extented to branches with
% SameLength=0, then the following lines needs to be changed (it's
% assumed that all branches have same length).
  schnitzNrs(i,:) = branches(i).schnitzNrs;
end
parent_schnitzes = schnitzNrs(:,1)';
%}

parent_schnitzes = ...
    arrayfun(@(branchIdx) branches(branchIdx).schnitzNrs(1), 1:numel(branches));

unique_parent_schnitzes = unique(parent_schnitzes);
disp(['There are ' num2str(length(unique_parent_schnitzes)) ' unique branches in data (parent cells: ' num2str(unique_parent_schnitzes) ')']);

%% DIVIDE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branch_groups = struct;
for group = 1:length(unique_parent_schnitzes)
  parent_cell = unique_parent_schnitzes(group);
  idx = find(parent_schnitzes == parent_cell);

  branch_groups(group).branches = branches(idx);
  branch_groups(group).parent_cell = parent_cell;
  branch_groups(group).nr_branches = length(idx);
end
