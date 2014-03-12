function [variation]=JK_variantie(x);

sum = 0;
count = 0;
for i = 1:length(x)
    var = ((x(i)-mean(x))^2);
    sum = sum + var;
    count = count + 1;
end
variation = sum / (count-1);