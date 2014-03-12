
function [output] = JK_corr_pooled(data, method, population_mean, population_stdev)
% 
% function [data,data_struct,data_array] = DJK_corr(data)
%
% Calculates the autocorrelation of some data
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
%
%

% if no method indicated, exit
if (~strcmp(method,'biased') & ~strcmp(method,'unbiased'))
    error('No method idicated');
end

% if too little data, exit
if (length(data) < 2), error('Not enough data'); end

% calc population_mean if not provided.
if exist('population_mean')~=1
    population_mean = 0;
    for n = 1:length(data) %Loop over each data point
        population_mean = population_mean + data(n);
    end
    population_mean = population_mean / length(data);
end

% calc of sum of all squares
sum_square = 0;
for n = 1:length(data) %Loop over each data point
    sum_square = sum_square + (data(n)-population_mean)*(data(n)-population_mean);
end

% if sum_square is zero, exit
if (sum_square == 0), error('sum_square is zero');, end

% calc of correlation
for delay = 0:length(data)-1 %Delay is 0 to length data - 1
    output(delay+1) = 0;
    for n = 1:length(data)-delay %Loop over each data point
        output(delay+1) = output(delay+1) + (data(n)-population_mean)*(data(n+delay)-population_mean);
    end
    output(delay+1) = output(delay+1)/sum_square;
end

% if unbiased
if (strcmp(method,'unbiased'))
    for delay = 0:length(data)-1
        output(delay+1) = output(delay+1) * (length(data) / (length(data) - delay) );
    end
end
