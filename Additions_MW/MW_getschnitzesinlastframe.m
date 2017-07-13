function lastSchnitzes = MW_getschnitzesinlastframe(schnitzcells)
% lastSchnitzes = MW_getlastschnitzes(schnitzcells)
%
% Function that returns the schnitzes that are present in the last frame, 
% based on frame_nrs field of the schnitzcells structure

% get last frame nr
lastFrame = max([schnitzcells(:).frame_nrs]);

% get schnitzes that exist in that last frame
lastSchnitzes = find(arrayfun(@(schnitzidx) any(schnitzcells(schnitzidx).frame_nrs==lastFrame), 1:numel(schnitzcells)))

end