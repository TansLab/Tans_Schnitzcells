% DJK_fillCellsToSize goes over all fields of structure, makes them cells
% if they are not yet, and enlarges cell to given cell_length
% * when enlarging the last value of cell is copied
%   so {[1 2] [3 4]} -> {[1 2] [3 4] [3 4] [3 4]}
% * if a field is a struct, it is ignored
%
% Code written by Daan Kiviet

function structure = DJK_fillCellsToSize(structure, cell_length);
% MAKE SURE THAT ALL FIELDS ARE CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structure = DJK_makeFieldsCells(structure);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ENLARGE CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldNames = fieldnames(structure);
for str = fieldNames'

  if isstruct(structure.(char(str))), continue; end
  
  oldCellSize = length(structure.(char(str)));
  if oldCellSize == 0
    warning(['Field ' char(str) ' is empty!']);
  else
    if cell_length > oldCellSize
      for i = oldCellSize+1:cell_length
        structure.(char(str)){i} = structure.(char(str)){oldCellSize};
      end
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
