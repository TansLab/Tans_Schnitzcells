function [power x0] = DJK_ExponentialFit( x, y);
% Performs an exponential fit by linearization of log values with base 2
%
%
p = polyfit(x,log2(y),1);
power = p(1);
x0 = 2^p(2);
