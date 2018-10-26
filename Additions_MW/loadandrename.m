
function output = loadandrename(theFilePath)
% Loads a parameter from .mat file, assumes there is only one parameter in
% there, and loads that (arbitrarily named) parameter into output.

% Load schnitzcells param
temp=load(theFilePath);
theParameterName{1} = fieldnames(temp);
theParameterNameString=theParameterName{1}{1};
output = temp.(theParameterNameString);
clear temp;
disp(['Loaded ' theParameterNameString ' into ''output'' parameter']);

end