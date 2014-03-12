function [meanMY]=getmeanMY(schnitzcells,schnitz)

MYdata = schnitzcells(schnitz).MY;
sum = 0;
count = 0;
for i=1:length(MYdata)
    temp = MYdata(i);
    if ~isnan( temp )
       sum = sum + temp;
       count = count + 1;
    end
end

%if count==0 return NaN;
meanMY = sum / count;

%beacause not all segmentation files contain fluorescence pictures. this
%poss a problem with other schnitzes, i should define which frames should
%be used. for instance: if frame number is n*11 then...)

%meanMY(i)=mean(MY)';
