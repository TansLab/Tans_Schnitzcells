 function [change_over_time] = JK_change_over_time(schnitzcells,startSchnitz,endSchnitz,datatype,population_mean);

    autocorrealtion = [];
 
    [data,data_struct,data_array] = DJK_getBranchData(schnitzcells, startSchnitz, endSchnitz); %to get data from specified branch
    
    [data_noNaN,data_struct_noNaN] = DJK_removeNaNData(data); % to delete schnitzes without fluorescence data
    
        data_tp = data_struct_noNaN.(datatype) - population_mean; 
           
        time = data_struct_noNaN.mins-data_struct_noNaN.mins(1);  
        
        change_over_time = [time;data_tp]';

end