


function L = MW_deletecellsattheedge(L)
% Cuts cells that are at the edge of the image. Every cell that is within
% MARGIN pixels of any edge, gets deleted.
%
% Input
%  - segmented image L
% Output
%  - segmented image L
% 

%% Parameters
MARGIN=20; % this should be equal to THICKNESS in MW_preprocessimagefadeedge

%% Function

if any(size(L) < MARGIN)
    error('L is too large for MW_deletecellsattheedge');
end

% Collect the values of the pixels that are at any edge
collection = [  L(1:end,MARGIN)',...
                L(1:end,end-MARGIN)',...
                L(MARGIN,1:end),...
                L(end-MARGIN,1:end) ]; 


% Now obtain the unique cell indices from that collection
cellIdxs = unique(collection);

% Determine where those cells are in the segmented image L
cellsToDelete = ismember(L,cellIdxs);

% Delete those cells 
L(cellsToDelete) = 0;

end





