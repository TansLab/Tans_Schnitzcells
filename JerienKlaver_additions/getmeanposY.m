function [meancenY]=getmeanposY(schnitzcells,schnitz)

cenYdata = schnitzcells(schnitz).cenY;
sum = 0;
count = 0;
for i=1:length(cenYdata)
    temp = cenYdata(i);
    if ~isnan( temp )
       sum = sum + temp;
       count = count + 1;
    end
end

%if count==0 return NaN;
meancenY = sum / count;

