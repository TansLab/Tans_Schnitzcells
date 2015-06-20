function outputString = escape(inputString)
% MW 2015/06
% Code stolen from: 
% http://stackoverflow.com/questions/15486611/how-to-add-before-all-special-characters-in-matlab

special = '[]{}()=''.().....,;:%%{%}!@_';

outputString = '';
for l = inputString
    if (length(find(special == l)) > 0)
        outputString = [outputString, '\', l];
    else
        outputString = [outputString, l];
    end
end

end