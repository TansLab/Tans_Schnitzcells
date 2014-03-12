
function [crosscorrelation_all] = JK_crosscor_all(schnitzcells,branchdata,method,datatype1,pop_mean1,pop_stdev1,datatype2,pop_mean2,pop_stdev2);

%crosscorrelation_all = struct();
%field_names = branchdata(1,:); %selects the upper row of the branchdata array


[crosscorrelation_all]=[];
for k = 1:length(branchdata(1,:)); %Loop over field names
    for i = branchdata(1,k);
        for j = branchdata(2,k);
            %disp( [i j]);
            temp = JK_getcrosscor(schnitzcells,i,j,method,datatype1,pop_mean1,pop_stdev1,datatype2,pop_mean2,pop_stdev2);
        end
        crosscorrelation_all{branchdata(1,k)} = temp;
    end
end