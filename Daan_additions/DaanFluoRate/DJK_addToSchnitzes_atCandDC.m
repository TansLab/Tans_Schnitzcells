% DJK_addToSchnitzes_atCandDC 
% same as DJK_addToSchnitzes_atCandDC but for C_time instead of C_time

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
function [schnitzcells] = DJK_addToSchnitzes_atCandDC(schnitzcells, field) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
field_atC                     = [field '_atC'];
field_atdC                    = [field '_atdC'];
% if field already existed, empty it
if isfield(schnitzcells,field_atC)
    schnitzcells=rmfield(schnitzcells,field_atC);
end
if isfield(schnitzcells,field_atdC)
    schnitzcells=rmfield(schnitzcells,field_atdC);
end
schnitzcells(1).(field_atC)   = [];
schnitzcells(1).(field_atdC)  = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);

  %------------------------------------------------------------------------
  % LOOP OVER ALL FLUOR MEASUREMENTS
  %------------------------------------------------------------------------
  for f = 1:length(s.C_time)
    age                   = find(s.time == s.C_time(f));
    s.(field_atC)(end+1)  = s.(field)(age);
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.C_time) > 1
    for f = 1:length(s.C_time)-1
      age                   = find(s.time == s.C_time(f));
      s.(field_atdC)(end+1) = s.(field)(age);
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if length(s.C_time)>0 % sometimes 1st schnitz does not have fluor images
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.C_time) > 0 & length(sE.C_time) > 0)
          age                   = find(s.time == s.C_time(end));
          s.(field_atdC)(end+1) = s.(field)(age);
        end
      end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
