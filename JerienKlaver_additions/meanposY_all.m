
function [meancenY_all] = meanposY_all(schnitzcells)

for i = 1:length(schnitzcells)
    tempList(i) = getmeanposY(schnitzcells,i);
  end
[meancenY_all] = [tempList]';

