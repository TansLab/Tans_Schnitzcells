% DJK_makeFieldsCells makes sure that all fields of given structure are cells
% * converts non-cells to cell: [1 2 3] -> {[1 2 3]}
% * fields that are structures are not converted to cells
%
% Code written by Daan Kiviet

function structure = DJK_makeFieldsCells(structure);

% LOOP OVER FIELDS AND IF NOT CELL ALREADY MAKE CELL %%%%%%%%%%%%%%%%%%%%%%
fieldNames = fieldnames(structure);
for str = fieldNames'
  if ~iscell(structure.(char(str))) && ~isstruct(structure.(char(str)))
    structure.(char(str)) = {structure.(char(str))};
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
