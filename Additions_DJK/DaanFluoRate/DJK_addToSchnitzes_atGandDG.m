% DJK_addToSchnitzes_atGandDG 
% same as DJK_addToSchnitzes_atYandDY but for G_time instead of Y_time
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
function [schnitzcells] = DJK_addToSchnitzes_atGandDG(schnitzcells, field) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
field_atG                     = [field '_atG'];
field_atdG                    = [field '_atdG'];
% if field already existed, empty it
if isfield(schnitzcells,field_atG)
    schnitzcells=rmfield(schnitzcells,field_atG);
end
if isfield(schnitzcells,field_atdG)
    schnitzcells=rmfield(schnitzcells,field_atdG);
end
schnitzcells(1).(field_atG)   = [];
schnitzcells(1).(field_atdG)  = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);

  %------------------------------------------------------------------------
  % LOOP OVER ALL FLUOR MEASUREMENTS
  %------------------------------------------------------------------------
  for f = 1:length(s.G_time)
    age                   = find(s.time == s.G_time(f));
    s.(field_atG)(end+1)  = s.(field)(age);
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.G_time) > 1
    for f = 1:length(s.G_time)-1
      age                   = find(s.time == s.G_time(f));
      s.(field_atdG)(end+1) = s.(field)(age);
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if length(s.G_time)>0 % sometimes 1st schnitz does not have fluor images
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.G_time) > 0 & length(sE.G_time) > 0)
          age                   = find(s.time == s.G_time(end));
          s.(field_atdG)(end+1) = s.(field)(age);
        end
      end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
