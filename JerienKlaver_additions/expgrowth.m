%determine exponential growth rate during whole cell cycle 
%28-11-2007 Jerien Klaver
function [mu] = expgrowth(schnitzcells,schnitz)
for i=1:length(schnitz)
    x=schnitzcells(1,schnitz(i)).mins; 
    y=schnitzcells(1,schnitz(i)).len;
    
    end
       
fit = polyfit(x, log(y), 1);
%to determine mu (which is the generation time per hour (mu=1 means a
%doubling time of 60 minutes):
mu(i)=[fit(1)*60/log(2)];

%to plot the linearized function use: 
%semilogy(x, y, 'o', x, exp(fit(2)).*exp(fit(1)*x))