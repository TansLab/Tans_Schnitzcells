function y = str3(x);
%
% function y = str3(x);
%
% A handy utility to convert an integer in the range 0-999 into a 
% 3-digit zero-padded string.  For example:
% 
%   str3(10) = '010'
%   str3(8)  = '008'
%  

y = num2str(x,'%5.3d');
