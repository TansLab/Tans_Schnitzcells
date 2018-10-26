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
% 'X'                 char coding for fluor color (e.g. G, Y, R, C)
%
% slighlty modified such that phaseFIelds are overwritten and not added one
% after the other if function is called multiple times (NW 2012-09-11)
%
function [schnitzcells] = MW_addToSchnitzes_atXandDX(schnitzcells, field,X) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
field_atX                     = [field '_at' X];
field_atdX                    = [field '_atd' X];
% if field already existed, empty it
if isfield(schnitzcells,field_atX)
    schnitzcells=rmfield(schnitzcells,field_atX);
end
if isfield(schnitzcells,field_atdX)
    schnitzcells=rmfield(schnitzcells,field_atdX);
end
schnitzcells(1).(field_atX)   = [];
schnitzcells(1).(field_atdX)  = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
timeField = [X '_time'];
for i = 1:length(schnitzcells)
  s = schnitzcells(i);

  %------------------------------------------------------------------------
  % LOOP OVER ALL FLUOR MEASUREMENTS
  %------------------------------------------------------------------------
  for f = 1:length(s.(timeField))
    age                   = find(s.time == s.(timeField)(f));
    s.(field_atX)(end+1)  = s.(field)(age);
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.(timeField)) > 1
    for f = 1:length(s.(timeField))-1
      age                   = find(s.time == s.(timeField)(f));
      s.(field_atdX)(end+1) = s.(field)(age);
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if length(s.(timeField))>0 % sometimes 1st schnitz does not have fluor images
      if (s.D ~= 0), sD = schnitzcells(s.D); end
      if (s.E ~= 0), sE = schnitzcells(s.E); end
      if (exist('sE') & exist('sD'))
        if (length(sD.(timeField)) > 0 & length(sE.(timeField)) > 0)
          age                   = find(s.time == s.(timeField)(end));
          s.(field_atdX)(end+1) = s.(field)(age);
        end
      end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
