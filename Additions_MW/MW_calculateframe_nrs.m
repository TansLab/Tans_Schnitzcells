function schnitzcells = MW_calculateframe_nrs(schnitzcells)
% 
% The Schnitzcells structure has an index for the frames with the images of
% the bacterial oclony. This index used to be called "frames", but this
% index was shifted +1. This has been removed, and to avoid conflicts, the
% index now has a new name, which is "frame_nrs". 
% Some old schnitzcells files still use the "frames" as index, and this
% function - very straightforwardly - converts the old index to the new,
% and adds it to the schnitzcells structure.
% See also: MW_addframenrsfix
% 
% INPUT
% Schnitzcells struct
%
% OUTPUT
% Schnitzcells struct

% backwards compatibility
if isfield(schnitzcells,'frames')
    
    for i = 1:numel(schnitzcells)           
        schnitzcells(i).frame_nrs = schnitzcells(i).frames-1; % N+1 bug    
    end
    
end

end