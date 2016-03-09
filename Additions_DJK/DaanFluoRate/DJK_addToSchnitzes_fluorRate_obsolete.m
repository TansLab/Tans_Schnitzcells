% DJK_addToSchnitzes_fluorRate takes a schnitzcells data set containing
% fluor data, calculates the fluor rate for all schnitzes and returns the
% schnitzcells data set with the rate.
%
% Sanders suggestion: STORES AT Y_TIME !!!!
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
%
function [schnitzcells] = DJK_addToSchnitzes_fluorRate(schnitzcells) 


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
schnitzcells(1).dY5_sum         = [];
schnitzcells(1).dY5_sum_dt      = [];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);
  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.Y_time) > 1
    for f = 1:length(s.Y_time)-1
      age                 = find(s.time == s.Y_time(f));
      age2                = find(s.time == s.Y_time(f+1));

      s.dY5_sum(end+1)    = s.Y5_sum_all(age2) - s.Y5_sum_all(age);
      s.dY5_sum_dt(end+1) = (s.Y5_sum_all(age2) - s.Y5_sum_all(age)) / (s.Y_time(f+1) - s.Y_time(f));
    end
  else
    f = 0;
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FOR LAST FLUOR MEASUREMENT CHECK WHETHER CHILDREN EXIST AND HAVE FLUOR
  %------------------------------------------------------------------------
  clear sD sE;
  if (s.D ~= 0), sD = schnitzcells(s.D); end
  if (s.E ~= 0), sE = schnitzcells(s.E); end
  if (exist('sE') & exist('sD'))
    if (length(sD.Y_time) > 0 & length(sE.Y_time) > 0)
      age                 = find(s.time == s.Y_time(end));
      age2                = find(sD.time == sD.Y_time(1));
      s.dY5_sum(end+1)    = sD.Y5_sum_all(age2) + sE.Y5_sum_all(age2) - s.Y5_sum_all(age);
      s.dY5_sum_dt(end+1) = (sD.Y5_sum_all(age2) + sE.Y5_sum_all(age2) - s.Y5_sum_all(age)) / (sD.Y_time(1) - s.Y_time(f+1));
    end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
%--------------------------------------------------------------------------


