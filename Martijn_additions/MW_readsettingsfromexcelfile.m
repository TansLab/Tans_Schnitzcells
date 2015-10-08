function settings = MW_readsettingsfromexcelfile(configfilepath)
    % read in excel file
    [ndata, text, alldata] = xlsread(configfilepath,'Configuration','A14:B45');

    % Create the variables that are listed in the configuration file
    % ===
    % alldata{i, 1} contains the desired names of the parameters, alldata{i, 2}
    % contains the values of the parameter.
    for i = 1:size(alldata,1)
        if ~(isempty(alldata{i, 1}) || isempty(alldata{i, 2}) || any(isnan(alldata{i, 2})) )
            % Tell user
            disp(['Processing ' alldata{i, 1}]);

            % create a parameter with the name contained by alldata{i, 1} 
            % (= left column of excel sheet), with value parameterValue.                
            % ===
            parameterValue = alldata{i, 2}; % put value in a var
            if isnumeric(parameterValue) 
                % if it is a numeric value, convert it to string w. numeric
                % value
                disp('Numeric');
                parameterValue = mat2str(parameterValue);
            elseif (~isempty(regexpi(alldata{i, 2}, '^[0123456789.]*$')))
                % If it is a string with a number, leave unchanged
                disp('Number');
            elseif alldata{i, 2}(1) == '['
                % If it is a string with a vector, leave unchanged
                disp('Vector');           
            else
                % If it is a string, add escape quotes such that it remains a
                % string
                disp('String');
                parameterValue = ['''' parameterValue ''''];
            end
            % Use the eval command to produce a parameter with the current value
            command = ['settings.' alldata{i, 1} '=' parameterValue];
            eval (command); 
        end
    end    

    disp('All configuration settings stored in ''settings'' struct.');
end