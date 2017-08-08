
function schnitzcells = MW_addToSchnitzes_fieldYatTime(schnitzcellsPath,timeFieldSelection, timeFieldSource, sourceField, someSuffix)%,saveTheFile)


% For a schnitzcells struct:
% This function creates a field for vectors with selected values of
% parameter Y, where selection is based on time vectors, timeFieldSelection
% is a vector with elements that are a subset of timeFieldSource. The
% resulting field will be named sourceField_someSuffix.

load(schnitzcellsPath);

%%
for schnitzIdx = 1:numel(schnitzcells)
    %%
    % Find out where we can find t1, t2, etc from the selection field
    % (eg production) in the source field (eg concentration)
    [sharedVals,idxsIntoA] = intersect(schnitzcells(schnitzIdx).(timeFieldSource),schnitzcells(schnitzIdx).(timeFieldSelection),'stable');

    schnitzcells(schnitzIdx).([sourceField someSuffix]) = ... 
        schnitzcells(schnitzIdx).(sourceField)(idxsIntoA');
    
    if numel(schnitzcells(schnitzIdx).(timeFieldSelection)) ~= numel(schnitzcells(schnitzIdx).([sourceField someSuffix]))
        error(['Something went wrong at schnitz #' num2str(schnitzIdx)]);
    end
end

if numel([schnitzcells(:).(timeFieldSelection)]) ~= numel([schnitzcells(:).([sourceField someSuffix])])
    error('Something went horribly wrong..');
end

%%
%if saveTheFile
save(schnitzcellsPath, 'schnitzcells');
%end

disp(['Done creating field ' sourceField someSuffix '.']);