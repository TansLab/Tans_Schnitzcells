
function schnitzcells = MW_addToSchnitzes_indices_atXdX(schnitzcellsPath,fluorColor)%,saveTheFile)
% function MW_addToSchnitzes_indices_atXdX(schnitzcellsPath,fluorColor)
% 
% Function loads and saves a schnitzcells parameter from schnitzcellsPath
% and creates fields that look like [0 0 0 1 0 0 0 1 0 0 0 1 ...] where
% ones indicate frames where fluorescence concentration (indices_at_X) or
% fluorescence production (indices_at_dX) where determined. The X is
% replaced by input parameter fluorColor.
%
% Input parameters:
%   > schnitzcellsPath
%   > fluorColor
%   > saveTheFile   
%
% Output parameters:
%   > schnitzcells - (note that schnitzcells struct is also saved.) with
%                     additional fields 
%           - schnitzcells.indices_at_X and 
%           - schnitzcells.indices_at_dX.

load(schnitzcellsPath);

for schnitzIdx = 1:numel(schnitzcells)
   schnitzcells(schnitzIdx).(['indices_at_'  fluorColor]) = ... 
        ismember(schnitzcells(schnitzIdx).time,schnitzcells(schnitzIdx).(['time_at' fluorColor]));
    
   schnitzcells(schnitzIdx).(['indices_at_d' fluorColor]) = ... 
        ismember(schnitzcells(schnitzIdx).time,schnitzcells(schnitzIdx).(['time_atd' fluorColor]));
    
   if isfield(schnitzcells,['d' fluorColor '5_time'])
       schnitzcells(schnitzIdx).(['indices_at_d' fluorColor '_time']) = ... 
            ismember(schnitzcells(schnitzIdx).time,schnitzcells(schnitzIdx).(['d' fluorColor '5_time']));
   end

end

%if saveTheFile
save(schnitzcellsPath, 'schnitzcells');
%end

disp(['Done adding fields indices_at_' fluorColor ' and indices_at_d' fluorColor '.']);