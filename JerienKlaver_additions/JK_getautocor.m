 function [autocorrelation] = JK_getautocor(schnitzcells,startSchnitz,endSchnitz,method,datatype,population_mean);

    autocorrealtion = [];
 
    [data,data_struct,data_array] = DJK_getBranchData(schnitzcells, startSchnitz, endSchnitz); %to get data from specified branch
    
    [data_noNaN,data_struct_noNaN] = DJK_removeNaNData(data); % to delete schnitzes without fluorescence data
    
        cov = DJK_corr((data_struct_noNaN.(datatype)),method,population_mean); %to determine covariance of fluorescence data in specified branch
    
         %autoc = cov((round(length(cov)/2):end))/(cov(round(length(cov)))/2); %to select only the second half of covairiance matrix 
         %(we are not intereseted in relation between fluorescence of t=0 and t= (0-x)
         %and to normalize values to 1
        autocorrelation = [(data_struct_noNaN.mins-data_struct_noNaN.mins(1))',cov'];  %to get matrix with time and autocorrelation

end