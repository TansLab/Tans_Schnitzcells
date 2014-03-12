function [N] = getcellsperschnitz(schnitzcells)
% 
% 
% Jerien Klaver 24 Oct 2007
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------
if (nargin < 1) 
  errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> convert_schnitzcells:',...
      '    Invalid input arguments.',...
      '    Try "help convert_schnitzcells".');
  error(errorMessage);
end

for i = 1:length(schnitzcells) %Loop over each old schnitz
    tempList(i) = schnitzcells(i).N;
end

% Save new list to memory
N = tempList;

return;