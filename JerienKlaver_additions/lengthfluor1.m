function [growth] = lengthfluor(schnitzcells,me)
% GETMY   ....
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


% frames -1, makes the frame number comparable to segmentation number.
    growth=[schnitzcells(me).frames-1;schnitzcells(me).len;..;100*schnitzcells(me).MY]'
