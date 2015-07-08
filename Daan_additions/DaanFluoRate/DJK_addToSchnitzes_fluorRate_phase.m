function [schnitzcells] = DJK_addToSchnitzes_fluorRate_phase(p,colorfield,myindex) 
% function [schnitzcells] = DJK_addToSchnitzes_fluorRate_phase(p,colorfield,myindex) 
% generalized to any input color by NW 2012-04 
%
% DJK_addToSchnitzes_fluorRate takes a schnitzcells data set containing
% fluor data, calculates the fluor rate for all schnitzes and returns the
% schnitzcells data set with the rate.
% The rate is determined by taking into account two points, and generally
% follows the formula: dY(t+dt)-dY(t). Additional fields are divided by dt
% and/or normalized. Time of reference t corresponds to the 'X_time' field,
% where X is the fluor (e.g. G, Y, C, R, ..).
% Note that there are two functions that determine the rate, this one and 
% PN_fluorRate_X. The latter has a different calculation, and takes into
% account t-1, t, t+1, thus making the time of reference of the two
% functions different. PN_fluorRate_X also adds the field dX5_time (this 
% fn doesn't do that) from which it can be seen which function was used.
% (-MW '15/04)

%
% IMPORTANT NOTE
% It is important to note that the rate determined here, is the absolute
% change in fluorescence per unit time per cell, i.e. NOT per volume 
% (or surface). Theoretically, it is expected to be always positive.
%
% Sanders suggestion: STORES AT Y_TIME !!!!
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'lengthField'       name of length Field
% 'myindex'           string allowing you to give a custom name to this field

if ~ischar(myindex)
    error('myindex input parameter should be string');
end

schnitzname = [p.tracksDir,p.movieName,'-Schnitz.mat'];
load(schnitzname);


% generate field names
% For different ways of determining fluor; colorfield tells you which fluor
% was used, and myindex allows manually determined different fields.
d_sum = ['d' upper(colorfield) myindex '_sum'];
d_sum_dt = ['d' upper(colorfield)  myindex '_sum_dt'];
d_sum_dt_ph = ['d' upper(colorfield)  myindex '_sum_dt_ph'];
d_sum_dt_vol = ['d' upper(colorfield)  myindex '_sum_dt_vol'];
sum_all = [upper(colorfield)  myindex '_sum_all'];

%d_sum = ['d' upper(colorfield) '5_sum'];
%d_sum_dt = ['d' upper(colorfield) '5_sum_dt'];
%d_sum_dt_ph = ['d' upper(colorfield) '5_sum_dt_ph'];
%d_sum_dt_vol = ['d' upper(colorfield) '5_sum_dt_vol'];
%sum_all = [upper(colorfield) '5_sum_all'];

fluortime = [upper(colorfield) '_time'];
phase_atfluor= ['phase_at' upper(colorfield)]; % phase = where in cell cycle 

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
ls = length(schnitzcells);
for ii = 1:ls
schnitzcells(ii).(d_sum)         = [];
schnitzcells(ii).(d_sum_dt)      = [];
schnitzcells(ii).(d_sum_dt_ph)   = [];
schnitzcells(ii).(d_sum_dt_vol)  = [];

schnitzcells(ii).(phase_atfluor)       = []; 
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);
  %------------------------------------------------------------------------
  % LOOP OVER FLUOR MEASUREMENTS EXCEPT LAST
  %------------------------------------------------------------------------
  if length(s.(fluortime))==0 % e.g. weird long cells which divide very fast 
                              % have sometimes no fluor image
      continue  % (NW 2012-08)
  else
  
      if length(s.(fluortime)) > 1
        % Loop over measured fluor values 
        for f = 1:length(s.(fluortime))-1
          % Note that age is normalized 
          age                 = find(s.time == s.(fluortime)(f)); % Age now
          age2                = find(s.time == s.(fluortime)(f+1)); % Age at next point

          % different ways of determining rate
          % absolute point to point difference
          s.(d_sum)(end+1)         =  s.(sum_all)(age2) - s.(sum_all)(age);
          % dY/dt
          s.(d_sum_dt)(end+1)      = (s.(sum_all)(age2) - s.(sum_all)(age)) /  (s.(fluortime)(f+1) - s.(fluortime)(f));
          % phase corrected
          s.(d_sum_dt_ph)(end+1)   = (s.(sum_all)(age2) - s.(sum_all)(age)) / ((s.(fluortime)(f+1) - s.(fluortime)(f))*(1 + s.phase(age)));
          % per volume
          s.(d_sum_dt_vol)(end+1)  = (s.(sum_all)(age2) - s.(sum_all)(age)) / ((s.(fluortime)(f+1) - s.(fluortime)(f))*(s.rp_volume(age)));

          s.(phase_atfluor)(end+1)       =  s.phase(age);
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
        if (length(sD.(fluortime)) > 0 & length(sE.(fluortime)) > 0)
         if i==1 % enable e.g. if first schnitzcell has no fluor at all

          s.(d_sum)(end+1)         =  0;
          s.(d_sum_dt)(end+1)      =  0;
          s.(d_sum_dt_ph)(end+1)   = 0 ;
          s.(d_sum_dt_vol)(end+1)  = 0;

          s.(phase_atfluor)(end+1)       =  0;
         else
          age                 = find(s.time == s.(fluortime)(end));
          age2                = find(sD.time == sD.(fluortime)(1));

          s.(d_sum)(end+1)         =  sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age);
          s.(d_sum_dt)(end+1)      = (sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age)) / (sD.(fluortime)(1) - s.(fluortime)(f+1));
          s.(d_sum_dt_ph)(end+1)   = (sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age)) / ((sD.(fluortime)(1) - s.(fluortime)(f+1))*(1 + s.phase(age)));
          s.(d_sum_dt_vol)(end+1)  = (sD.(sum_all)(age2) + sE.(sum_all)(age2) - s.(sum_all)(age)) / ((sD.(fluortime)(1) - s.(fluortime)(f+1))*(s.rp_volume(age)));

          s.(phase_atfluor)(end+1)       =  s.phase(age);
         end
        end
      end
      %------------------------------------------------------------------------

      schnitzcells(i) = s;
  end % exclude weird cells
end 

save(schnitzname,'schnitzcells');