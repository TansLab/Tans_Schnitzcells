function [cenX] = getabsX(schnitzcells)
% GET absolute value of X coordinate   ....
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
    tempList(i) = schnitzcells(i).cenX;
end

% Save new list to memory
cenX = tempList;

return;