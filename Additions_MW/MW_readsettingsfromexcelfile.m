function [settings, alldata] = MW_readsettingsfromexcelfile(settings)
    % Input argument required:
    % struct settings that has a field called configfilepath, which holds
    % the location of the configfile.


    % configfilepath = CONFIGFILEPATH

    %% read in excel file (note that we start reading at line 14)
    [ndata, text, alldata] = xlsread(settings.configfilepath,'Configuration','A14:B100');

    %% Create the variables that are listed in the configuration file
    % ===
    % alldata{i, 1} contains the desired names of the parameters, alldata{i, 2}
    % contains the values of the parameter.
    processedFieldNames = {};
    for i = 1:size(alldata,1)
        if ~(isempty(alldata{i, 1}) || isempty(alldata{i, 2}) || any(isnan(alldata{i, 2})) )
            
            % put value in a var
            parameterValue = alldata{i, 2}; 
            
            % Tell user
            disp(['Processing ' alldata{i, 1}]);
            disp(['Value = ' parameterValue ]);
            
            % Double check whether it's not double (leads to hairy situations)
            if any(strcmp(processedFieldNames,alldata{i, 1}))
                error('There are duplicate field names in Excel file, please correct!')
            end
            
            %% create a parameter with the name contained by alldata{i, 1} 
            % (= left column of excel sheet), with value parameterValue.                
            % ===            
            if isnumeric(parameterValue) 
                % if it is a numeric value, convert it to string w. numeric
                % value
                disp('Numeric');
                parameterValue = mat2str(parameterValue);
            elseif (~isempty(regexpi(alldata{i, 2}, '^[0123456789.]*$')))
                % If it is a string with a number, leave unchanged
                disp('Number');
            elseif alldata{i, 2}(1) == '[' || alldata{i, 2}(1) == '{'
                % If it is a string with a vector or cell, leave unchanged
                disp('Vector or cell');           
            else
                % If it is a string, add escape quotes such that it remains a
                % string
                disp('String');
                parameterValue = ['''' parameterValue ''''];
            end
            
            %% Use the eval command to produce a parameter with the current value
            command = ['settings.' alldata{i, 1} '=' parameterValue ';']
            eval (command); 
            
            % Administration
            processedFieldNames{end+1} = alldata{i, 1};
        end
    end    

    disp('All configuration settings stored in ''settings'' struct.');
end