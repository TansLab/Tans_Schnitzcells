% DJK_divide_branch_data divides branch data in seperate groups, each
% originating from the same parent cell. 
%
% returns a collection of branches (a struct of structs)
%
function branch_groups = DJK_divide_branch_data(branches);

% FIRST CHECK HOW MANY BRANCHES IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(branches)
  schnitzNrs(i,:) = branches(i).schnitzNrs;
end
parent_schnitzes = schnitzNrs(:,1)';
unique_parent_schnitzes = unique(parent_schnitzes);
disp(['There are ' num2str(length(unique_parent_schnitzes)) ' unique branches in data (parent cells: ' num2str(unique_parent_schnitzes) ')']);

% DIVIDE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branch_groups = struct;
for group = 1:length(unique_parent_schnitzes)
  parent_cell = unique_parent_schnitzes(group);
  idx = find(parent_schnitzes == parent_cell);

  branch_groups(group).branches = branches(idx);
  branch_groups(group).parent_cell = parent_cell;
  branch_groups(group).nr_branches = length(idx);
end
