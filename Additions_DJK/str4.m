function y = str4(x)
%
% function y = str4(x);
%
% A handy utility to convert an integer in the range 0-99 into a 
% 4-digit zero-padded string.  For example:
% 
%   str2(10) = '0010'
%   str2(8)  = '0008'
%  

y = num2str(x,'%5.4d');