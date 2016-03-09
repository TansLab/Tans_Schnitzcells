
% This functionality has been moved to NW_detecSlowSchnitzes.
% (File kept if necessary for other purposes.)














%%

associatedFieldNames =  {'G_time','G6_mean', 'muP11_fitNew'};

%% Find schnitzes that have NaN growth rates (probably negatively growing)
hits = [];
for i= 1:numel(schnitzcells)
    idxWithNans = find(isnan([schnitzcells(i).(associatedFieldNames{3})]));
    if ~isempty(idxWithNans)
        hits(end+1) = i;
    end
end

hits

%%

disp('Showing vectors per hit:');
disp('==================================================================');
schnitzcells(hits).muP11_fitNew
disp('==================================================================');
disp('Showing badSchnitzes:');
mat2str(hits)