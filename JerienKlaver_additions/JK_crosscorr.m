
function [output] = JK_crosscorr(data1, data2, method, pop_mean1,pop_stdev1,pop_mean2,pop_stdev2);
% 
% function [data,data_struct,data_array] = DJK_corr(data)
%
% Calculates the crosscorrelation of some data
% 
% method can be:
%   * 'biased'
%       The formula used for the autocorrelation function is essentially 
%       the same as the one is Nature439_608.
%   * 'unbiased'
%       The formula used for the autocorrelation function is essentially
%       the same as the one on the wikipedia page.
%
% As a normalisation, the population_mean will be substracted from the
% data. If the population_mean is not given, it will be calculated as the
% mean of your data set. NOTE, that if the mean of your data set is
% different from the real mean, the output will not give the correct
% correlation
% As a normalisation, the population_stdev will be substracted from the
% data. If the population_mean is not given, it will be calculated as the
% mean of your data set.
%

% if no method indicated, exit
if (~strcmp(method,'biased') & ~strcmp(method,'unbiased'))
    error('No method idicated');
end

% if too little data, exit
if (length(data1) < 3), error('Not enough data'); end

% calc population1_mean if not provided.
if exist('pop_mean1')~=1
    pop_mean = 0;
    for n = 1:length(data1) %Loop over each data point
        pop_mean = pop_mean + data1(n);
    end
    pop_mean1 = pop_mean / length(data1);
end

% calc population2_mean if not provided.
if exist('pop_mean2')~=1
    pop_mean = 0;
    for n = 1:length(data2) %Loop over each data point
        pop_mean = pop_mean + data2(n);
    end
    pop_mean2 = pop_mean / length(data2);
end

% calc pop_stdev1 if not provided.
if exist('pop_stdev1')~=1
    pop_SS = 0;
    for n = 1:length(data1) %Loop over each data point
        pop_SS = pop_SS + (data1(n)- pop_mean1)^2;
    end
    pop_stdev1 = sqrt(pop_SS / length(data1));
end

% calc pop_stdev2 if not provided.
if exist('pop_stdev2')~=1
    pop_SS = 0;
    for n = 1:length(data2) %Loop over each data point
        pop_SS = pop_SS + (data2(n)- pop_mean2)^2;
    end
    pop_stdev2 = sqrt(pop_SS / length(data2));
end

% calc of correlation 
% time-frame shifts
for delay = 0:length(data1)-1 %Delay is 0 to length data - 1
    output(delay+1) = 0;
    for n = 1:length(data1)-delay %Loop over each data point
        output(delay+1) = output(delay+1) + (data1(n)-pop_mean1)*(data2(n+delay)-pop_mean2);
    end
    output(delay+1) = output(delay+1)/(pop_stdev1*pop_stdev2*length(data1));
end

% if unbiased
if (strcmp(method,'unbiased'))
    for delay = 0:length(data1)-1
        output(delay+1) = output(delay+1) * length(data1) / (length(data1)-delay);
    end
end
end