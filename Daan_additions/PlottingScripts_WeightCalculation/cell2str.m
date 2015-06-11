function string = cell2str(cell)
if nargin~=1
    warning('CELL2STR:Nargin','Takes 1 input argument.');
end

if iscell(cell)
  string = '{ ';
  for i = 1:length(cell)
    string = [string num2str(cell{i}) ' '];
  end
  string = [string '}'];
else
  string = num2str(cell);
end
