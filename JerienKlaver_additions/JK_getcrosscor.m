function [crosscorrelation] = JK_getcrosscor(schnitzcells,startSchnitz,endSchnitz,method,datatype1,pop_mean1,pop_stdev1,datatype2,pop_mean2,pop_stdev2);

         
    [data,data_struct,data_array] = DJK_getBranchData(schnitzcells, startSchnitz, endSchnitz); %to get data from specified branch
    
    [data_noNaN,data_struct_noNaN] = DJK_removeNaNData(data); % to delete schnitzes without fluorescence data
    %to determine crosscorrelation of fluorescence data in specified branch
        
    cov = JK_crosscorr((data_struct_noNaN.(datatype1)),(data_struct_noNaN.(datatype2)),method,pop_mean1,pop_stdev1,pop_mean2,pop_stdev2); 
        %to elongate time array to the same size as cov array:
        %for i=1:length(data_struct_noNaN.mins)-1;
         %   A(i)=-data_struct_noNaN.mins(1+length(data_struct_noNaN.mins)-i);
        %end
        %for i=length(data_struct_noNaN.mins):length(data_struct_noNaN.mins)*2-1;
         %   A(i)= data_struct_noNaN.mins(i-length(data_struct_noNaN.mins)+1);
        %end
        A = data_struct_noNaN.mins-data_struct_noNaN.mins(1);
        crosscorrelation = [A',cov'];  %to get matrix with time and crosscorrelation

end