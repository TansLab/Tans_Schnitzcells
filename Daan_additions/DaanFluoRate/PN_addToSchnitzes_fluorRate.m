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
function [schnitzcells] = PN_addToSchnitzes_fluorRate(p,colorfield) 

schnitzname = [p.tracksDir,p.movieName,'-Schnitz.mat'];
load(schnitzname);

d_sum = ['d' upper(colorfield) '5_sum'];
d_sum_dt = ['d' upper(colorfield) '5_sum_dt'];
time = [upper(colorfield) '_time'];
sum_all = [upper(colorfield) '5_sum_all'];

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
schnitzcells(1).(d_sum)         = [];
schnitzcells(1).(d_sum_dt)      = [];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  
    schnitzcells(i).(d_sum)         = [];  %blubb %delete old data
    schnitzcells(i).(d_sum_dt)      = [];  %blubb
    
    s = schnitzcells(i);
  
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.(time)) > 1
    for f = 1:length(s.(time))-1
      age                 = find(s.time == s.(time)(f));
      age2                = find(s.time == s.(time)(f+1));

      s.(d_sum)(end+1)    = s.(sum_all)(age2) - s.(sum_all)(age);
      s.(d_sum_dt)(end+1) = (s.(sum_all)(age2) - s.(sum_all)(age)) / (s.(time)(f+1) - s.(time)(f));
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
    if (length(sD.(time)) > 0 & length(sE.(time)) > 0)
      age                 = find(s.time == s.(time)(end));
      age2                = find(sD.time == sD.(time)(1));
      s.(d_sum)(end+1)    = sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age);
      s.(d_sum_dt)(end+1) = (sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age)) / (sD.(time)(1) - s.(time)(f+1));
    end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 

save(schnitzname,'schnitzcells');

