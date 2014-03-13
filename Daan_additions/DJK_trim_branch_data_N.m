% DJK_trim_branch_data trims data from beginning of branches untill there
% are N main branches left.
%
function trimmed_branches = DJK_trim_branch_data_N(branches,N);

% GET BRANCH INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment NW 2013-04-11: if CrossCorrs will be extented to branches with
% SameLength=0, then the following two lines needs to be changed (it's
% assumed that all branches have same length).
for i=1:length(branches)
  schnitzNrs(i,:) = branches(i).schnitzNrs;
end

% LOOP OVER TIME POINTS AND DETERMINE NR OF BRANCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(schnitzNrs(1,:))
  parent_schnitzes = schnitzNrs(:,j)';
  unique_parent_schnitzes = unique(parent_schnitzes);
  disp(['At timepoint ' num2str(j) ' there are ' num2str(length(unique_parent_schnitzes)) ' unique branches in data (parent cells: ' num2str(unique_parent_schnitzes) ')']);
  if length(unique_parent_schnitzes) >= N
    break;
  end
end

% TRIM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trimmed_branches = branches;
fields = fieldnames(branches);
for f = 1:length(fields)
  for i=1:length(trimmed_branches)
    trimmed_branches(i).(char(fields(f))) = branches(i).(char(fields(f)))(j:end);
  end
end

