% DJK_getNrCellsAtEnd counts the number of cells in the last frame of the
% provided schnitzcells structure. Also returns the total length of these
% cells (so DJK_addToSchnitzes_length must have been run)
%
% OUTPUT
% 'count'           number of cells in the last frame
% 'total_length'    total length of these cells
%
% REQUIRED ARGUMENTS:
% 's'               schnitzcells structure
%  

function [count total_length] = DJK_getNrCellsAtEnd(s)

%--------------------------------------------------------------------------
% DO CALCULATION
%--------------------------------------------------------------------------
allFrames = unique([s(:).frames]);
lastFrame = max(allFrames);
count = 0;
total_length = 0;
for i = 1:length(s)
  idx = find(s(i).frames==lastFrame);
  if idx
    count = count + 1;
    total_length = total_length + s(i).length_fitNew(idx(1));
  end
end
disp(['At last frame ' str3(lastFrame) ' there were ' num2str(count) ' cells with a total length of ' num2str(round(total_length)) ' um.']);
%--------------------------------------------------------------------------
