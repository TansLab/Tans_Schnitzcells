
function schnitzcells = MW_addToSchnitzes_fieldYatX(schnitzcellsPath,fieldY, atX, someSuffix)%,saveTheFile)


warning('This function is untested!');
% It was not needed in the end because I used the timefield for selection
% in the problem where this function was required.



%
% For a schnitzcells struct:
% This function creates a field for vectors with selected values of
% parameter Y, where selection is based on binary vectors supplied by the
% field atX. The resulting field will be named fieldY_someSuffix.

load(schnitzcellsPath);

%%
for schnitzIdx = 1:numel(schnitzcells)
   schnitzcells(schnitzIdx).([fieldY someSuffix]) = schnitzcells(schnitzIdx).(fieldY)(schnitzcells(schnitzIdx).(atX));
end

%%
%if saveTheFile
save(schnitzcellsPath, 'schnitzcells');
%end

disp(['Done creating field ' fieldY someSuffix '.']);



