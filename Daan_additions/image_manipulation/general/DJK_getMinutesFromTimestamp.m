function time = DJK_getMinutesFromTimestamp(timestamp)
% function time = DJK_getMinutesFromTimestamp(timestamp);
%
% A utility to convert a timestamp into minutes
%  
[y,m,d,h,mi,s] = datevec(timestamp);    
time = s/60 + mi + h*60 + d*24*60; % + m*30*24*60;
