function [autocorrelation]=JK_autocorrelation(x);

sum = 0;
count = 0;
for i = 1:length(x)
    var = ((x(i)-mean(x))^2);
    sum = sum + var;
    count = count + 1;
end
variation = sum / count;

sum = 0;
count = 0;
a = length(x);
for i = 0:(length(x)-1)
    sum = 0;
    count = 0;
     for j = 1:a
     cor = (x(j)-mean(x))*(x(j+i)-mean(x));
     sum = sum + cor;
     count = count +1;
     end
    Ecorr = sum / count;
    autocorrelation(i+1) = Ecorr/variation;
    a = a-1;
end
autocorrelation(i+1) = Ecorr/variation;
end