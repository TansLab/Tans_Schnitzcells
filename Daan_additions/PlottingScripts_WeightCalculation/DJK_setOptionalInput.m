% DJK_setOptionalInput creates a field in given structure and set is to given value
% * example: DJK_setOptionalInput(fg, 'fg.lineXY.LineStyle', '--');
%
% Code written by Daan Kiviet

function structure = DJK_setOptionalInput(structure, input, value);

% CHECK WHETHER INPUT EXISTS AND IF NOT CREATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first extract original structure name from input, and create internally
idx = strfind(input,'.');
eval([input(1:idx(1)-1) ' = structure;']);

% test whether input exists now internally
try
  eval([input ';']);

% if input does not exist internally, create it, and set output to internal structure
catch
  eval([input ' = {value};']);
  eval(['structure = ' input(1:idx(1)-1) ';']);
  disp([' |---> set ' input ' to ' cell2str(value)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
