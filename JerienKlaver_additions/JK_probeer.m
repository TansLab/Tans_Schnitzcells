
%deze versie werkt nog niet. nu wil ik branchdata met meer dan 1 input en ik wil
%uiteindelijk een structure met per startschnitz een autocorrelatie. het
%gaat mis bij 20

%branchdata array should contain two rows with start (first row) end
%endschnitzes (second)

function [autocorrelation_all] = JK_probeer(schnitzcells,branchdata)

%autocorrelation_all = struct();
%field_names = branchdata(1,:); %selects the upper row of the branchdata array



for k = 1:length(branchdata(1,:)) %Loop over field names
    for i = branchdata(1,k)
        for j = branchdata(2,k)
            temp = JK_getautocor(schnitzcells,i,j);
        end
        %autocorrelation_all(k) = setfield(autocorrelation_all,'autocorrelation',temp);
        autocorrelation_all{branchdata(1,k)} = temp;
    end
end

