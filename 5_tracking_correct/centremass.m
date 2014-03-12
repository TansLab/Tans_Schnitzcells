function [x,y]= centremass(a)
% [x,y]= centremass(a)
%
% finds centre of mass of black and white image

a= +(a>0);

[fx, fy]= find(a>0);
tot= sum(sum(a));
x = sum(fx)/tot;
y = sum(fy)/tot;
