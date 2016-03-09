% DJK_addToSchnitzes_atRandDR 
% same as DJK_addToSchnitzes_atYandDY but for R_time instead of Y_time
% (->ribosome consturcts)
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'field'             name of field
%
% slighlty modified such that phaseFIelds are overwritten and not added one
% after the other if function is called multiple times (NW 2012-09-11)
%
function [schnitzcells] = DJK_addToSchnitzes_atRandDR(schnitzcells, field) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
field_atR                     = [field '_atR'];
field_atdR                    = [field '_atdR'];
% if field already existed, empty it
if isfield(schnitzcells,field_atR)
    schnitzcells=rmfield(schnitzcells,field_atR);
end
if isfield(schnitzcells,field_atdR)
    schnitzcells=rmfield(schnitzcells,field_atdR);
end
schnitzcells(1).(field_atR)   = [];
schnitzcells(1).(field_atdR)  = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);

  %------------------------------------------------------------------------
  % LOOP OVER ALL FLUOR MEASUREMENTS
  %------------------------------------------------------------------------
  for f = 1:length(s.R_time)
    age                   = find(s.time == s.R_time(f));
    s.(field_atR)(end+1)  = s.(field)(age);
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.R_time) > 1
    for f = 1:length(s.R_time)-1
      age                   = find(s.time == s.R_time(f));
      s.(field_atdR)(end+1) = s.(field)(age);
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if length(s.R_time)>0 % sometimes 1st schnitz does not have fluor images
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.R_time) > 0 & length(sE.R_time) > 0)
          age                   = find(s.time == s.R_time(end));
          s.(field_atdR)(end+1) = s.(field)(age);
        end
      end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
