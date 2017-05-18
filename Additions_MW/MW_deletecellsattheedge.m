


function L = MW_deletecellsattheedge(p,L)
% Cuts cells that are at the edge of the image. Every cell that is within
% MARGIN pixels of any edge, gets deleted.
%
% Input
%  - segmented image L
% Output
%  - segmented image L
% 

%% Parameters
MARGIN=30; % this should be tuned with size of filter in MW_preprocessimagefadeedge

%% Function

if any(size(L) < MARGIN)
    error('L is too large for MW_deletecellsattheedge');
end

% Collect the values of the pixels that are at any edge
collection1=[];collection2=[];collection3=[];collection4=[];
% depending on where the cells disappear
if any(p.mothermachine==3)
    collection1=L(1:end,1:MARGIN);          % left
end
if any(p.mothermachine==4)
    collection2=L(1:end,end-MARGIN:end);    % right
end
if any(p.mothermachine==2)
    collection3=L(1:MARGIN,1:end);          % top
end
if any(p.mothermachine==1)
    collection4=L(end-MARGIN:end,1:end);    % bottom
end

% Now obtain the unique cell indices from that collection
cellIdxs = unique( [collection1(:)' collection2(:)' collection3(:)' collection4(:)'] );

% Determine where those cells are in the segmented image L
cellsToDelete = ismember(L,cellIdxs);

% Delete those cells 
L(cellsToDelete) = 0;

% figure; PN_imshowlabel(p,L,[],[],[]);

%{
% show edges 
if any(p.mothermachine==3)
    L(1:end,1:MARGIN) =1;          % left
end
if any(p.mothermachine==4)
    L(1:end,end-MARGIN:end) =1;    % right
end
if any(p.mothermachine==2)    
    L(1:MARGIN,1:end) =1;          % top
end
if any(p.mothermachine==1)
    L(end-MARGIN:end,1:end) =1;    % bottom
end

PN_imshowlabel(p,L,[],[],[]);
%}

end





