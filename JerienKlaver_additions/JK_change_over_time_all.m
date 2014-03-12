
function [change_over_time_all] = JK_change_over_time_all(schnitzcells,branchdata,datatype,population_mean)
% branchdata is a two row array with in the first row the startnumber of a
% branch and in the second row the end-schnitz
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


for k = 1:length(branchdata(1,:)) %Loop over field names
    for i = branchdata(1,k)
        for j = branchdata(2,k)
            
            temp = JK_change_over_time(schnitzcells,i,j,datatype,population_mean);
        end
           change_over_time_all{branchdata(2,k)} = temp;
    end
end