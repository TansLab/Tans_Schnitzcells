function branches_nonempty = MW_remove_empty_branches(branches,N)
% uses the field schnitzNrs to determine which branches are empty and
% remove those
%
% arguments:
%   branches:   branch structure 
% optional arguments:
%   N:          minimum number of point in the branch that we keep

% 
if ~exist('N', 'var')
    N=1;
end

%% Count branches
numelBranches = arrayfun(@(x) numel(x.schnitzNrs),branches);

%% Select branches that are longer than N
branches_nonempty = branches(numelBranches>N);


%%
% Less elegant way:
%{
%%
branches_nonempty=struct;
allFields =  fieldnames(branches);
for fieldIdx = 1:numel(allFields)
    fieldName = allFields{fieldIdx};
    branches_nonempty.(fieldName) = [];
end

%%
count=0;
for branchIdx = 1:numel(branches)
    
   if numel(branches(branchIdx).schnitzNrs)>=N
       
       count=count+1;
       branches_nonempty(count) = branches(branchIdx);
       
   end

end
%}

end


