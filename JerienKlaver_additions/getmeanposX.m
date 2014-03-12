function [meancenX]=getmeanposX(schnitzcells,schnitz)

cenXdata = schnitzcells(schnitz).cenX;
sum = 0;
count = 0;
for i=1:length(cenXdata)
    temp = cenXdata(i);
    if ~isnan( temp )
       sum = sum + temp;
       count = count + 1;
    end
end

%if count==0 return NaN;
meancenX = sum / count;

