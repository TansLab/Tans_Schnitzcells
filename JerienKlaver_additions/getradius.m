function [r] = getradius(schnitzcells)
% GET X and Y coordinate 
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
% find position X and its deviation from the center (center is thought to
% be 750) and multiply by position Y and its deviation from the center.
for i = 1:length(schnitzcells) %Loop over each old schnitz
    tempList(i) = sqrt(((schnitzcells(i).cenX-750)^2)+(schnitzcells(i).cenY-750)^2);
end

% Save new list to memory
r = tempList;



return;