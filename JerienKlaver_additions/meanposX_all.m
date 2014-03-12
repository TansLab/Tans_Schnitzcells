
function [meancenX_all] = meanposX_all(schnitzcells)

for i = 1:length(schnitzcells)
    tempList(i) = getmeanposX(schnitzcells,i);
  end
[meancenX_all] = [tempList]';

