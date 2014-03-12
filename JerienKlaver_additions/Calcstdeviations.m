function [stdeviations] = Calcstdeviations(data,averages)


groupTime = round(data(1,1));
groupNumber = 1;
sum1 = 0;
sum2 = 0;
count = 0;

for i=1:length(data(:,1))
    if ( round(data(i,1)) ~= groupTime)
        stdeviations(groupNumber,2) = sqrt(sum1 / count);
        stdeviations(groupNumber,1)= sum2 / count;
        sum1 = 0;
        sum2 = 0;
        count = 0;
        groupTime = round(data(i,1));
        groupNumber = groupNumber + 1; 
    end
   
    sum1 = sum1 + (data(i,2)-averages(groupNumber,2))^2;
    sum2 = sum2 + round(data(i,1));
    count = count+1;
   
end
stdeviations(groupNumber,2) = sqrt(sum1 / count);
        stdeviations(groupNumber,1)= sum2 / count;
