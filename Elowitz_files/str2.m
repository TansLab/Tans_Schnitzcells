function y = str2(x);
%
% function y = str2(x);
%
% A handy utility to convert an integer in the range 0-99 into a 
% 1-digit zero-padded string.  For example:
% 
%   str2(10) = '10'
%   str2(8)  = '08'
%  

y = num2str(x,'%5.2d');
