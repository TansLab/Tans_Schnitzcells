% DJK_addToSchnitzes_atYandDY 
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
function [schnitzcells] = DJK_addToSchnitzes_atYandDY(schnitzcells, field) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
field_atY                     = [field '_atY'];
field_atdY                    = [field '_atdY'];
% if field already existed, empty it
if isfield(schnitzcells,field_atY)
    schnitzcells=rmfield(schnitzcells,field_atY);
end
if isfield(schnitzcells,field_atdY)
    schnitzcells=rmfield(schnitzcells,field_atdY);
end
schnitzcells(1).(field_atY)   = [];
schnitzcells(1).(field_atdY)  = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);

  %------------------------------------------------------------------------
  % LOOP OVER ALL FLUOR MEASUREMENTS
  %------------------------------------------------------------------------
  for f = 1:length(s.Y_time)
    age                   = find(s.time == s.Y_time(f));
    s.(field_atY)(end+1)  = s.(field)(age);
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.Y_time) > 1
    for f = 1:length(s.Y_time)-1
      age                   = find(s.time == s.Y_time(f));
      s.(field_atdY)(end+1) = s.(field)(age);
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if length(s.Y_time)>0 % sometimes 1st schnitz does not have fluor images
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.Y_time) > 0 & length(sE.Y_time) > 0)
          age                   = find(s.time == s.Y_time(end));
          s.(field_atdY)(end+1) = s.(field)(age);
        end
      end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
