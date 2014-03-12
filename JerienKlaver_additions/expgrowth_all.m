%determine exponential growth rate during whole cell cycle for all cells
%28-11-2007 Jerien Klaver
function [mu_all] = expgrowth_all(schnitzcells)

for i = 1:length(schnitzcells)
    tempList(i) = expgrowth(schnitzcells,i);
  end
[mu_all] = [tempList]';

