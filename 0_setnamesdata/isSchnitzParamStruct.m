function result = isSchnitzParamStruct(p)
%ISSHNITZPARAMSTRUCT Test if variable is a schnitzcells structure.
%
%   isShnitzParamStruct(p) returns the logical true (1) if p is a 
%   schnitzcells parameter structure, false (0) otherwise.
%

% Note: I wanted to use exist('p.movieName') but exist won't work on fields.

result = false;

if isstruct(p) 
  try
    name = p.movieName;
    result = true;
  catch
    % disp('p is missing field "movieName"');
  end
end
