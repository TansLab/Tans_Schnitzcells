%determine mean fluor during whole cell cycle for all cells
%28-11-2007 Jerien Klaver
function [meanMY_all] = meanfluor_all(schnitzcells)

for i = 1:length(schnitzcells)
    tempList(i) = getmeanMY(schnitzcells,i);
  end
[meanMY_all] = [tempList]';

