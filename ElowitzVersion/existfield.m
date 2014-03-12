function result = existfield(structure, fieldnameString)
%fieldexist Test if field in a structure exists
%
%   existfield(structure, fieldnameString) returns the logical true (1) if 
%   field (field given as a string) exist a field of the given structure, 
%   and false (0) otherwise.
%

% Note: I wanted to use exist('p.fieldname') but exist won't work on fields.

result = false;

if ~isstruct(structure)
  % structure doesn't exist, so return false
  disp('structure does not exist');
  result = false
  return
end

try
  eval(['testname = structure.' fieldnameString ';']);
  % disp(['field exists and has value ' testname]);
  result = true;
catch
  % disp(['structure field missing field ' fieldnameString]);
  result = false;
end

return 