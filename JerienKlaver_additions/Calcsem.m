function [sem] = Calcsem(data,averages,binwidth)
%use this function to calculate sSEM values of X in
%matrix with time in first collumn and X data in second. Depending on binwidth, the average is
%calculated over a certain bin of times. 

    %for exapmle: binwidth = 10
    %timebins from 0-5, 5-15 and so on
    
    
data(:,1) = data(:,1)/binwidth;
groupTime = binwidth*round(data(1,1));
groupNumber = 1;
sum1 = 0;
sum2 = 0;
count = 0;

for i=1:length(data(:,1))
    if (binwidth*round(data(i,1)) ~= groupTime)
        sem(groupNumber,2) = sqrt(sum1 / (count-1))/(sqrt(count));
        sem(groupNumber,1)= sum2*binwidth / count;
        sum1 = 0;
        sum2 = 0;
        count = 0;
        groupTime = binwidth*round(data(i,1));
        groupNumber = groupNumber + 1; 
    end
   
    sum1 = sum1 + (data(i,2)-averages(groupNumber,2))^2;
    sum2 = sum2 + round(data(i,1));
    count = count+1;
   
end
        sem(groupNumber,2) = sqrt(sum1 / count)/(sqrt(count));
        sem(groupNumber,1)= sum2*binwidth / count;
