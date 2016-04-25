function [cellid PP DD EE]= cellids(P,D,E,G);
% 
% function [cellid PP DD EE]= cellids(P,D,E,G);
%  used in recalc_schnitz()
% 
% calculates the label_id -> schnitznumber correspondence
% 
% Output
%   - cellid        cellid{framenr}(cellno) gives the corresponding lineage
%                   number (schnitznumber) for a given cell number (i.e.
%                   the number that cell has in the frame) for frame
%                   framenr.
%   - PP            PP(schnitzindx) gives the schnitzindex of the parent 
% 
% Input
%   - P     P{framenr}(i) holds cellnr of cell i in frame framenr, if it
%           did not divide. Otherwise 0.
%   - D     D{framenr}(i) holds cellnr of cell i in frame framenr, if it
%           did divide. Corresponding P2 holds cellno in next frame of 
%           daughter 1.
%   - E     E{framenr}(i) holds cellnr of cell i in frame framenr, if it
%           did divide. Corresponding P2 holds cellno in next frame of 
%           daughter 2.
%   - G     Ghost cells that have no parent?
% and parent (PP), daughter one (DD), and daughter two (EE)arrays
% labeled: <cell id>= <??>(cell_id)
% T.B. 04/05 

% MW:
% See also DJK_data_treat
% This function is used also in MW_makeSchnitzFileFromTracking, when 
% recalc_schnitz is called. 
% In MW_makeSchnitzFileFromTracking, 
% [P, D, E, G] = DJK_data_treat(p)
% , is executed to load the values. Note that these are all cells now, so
% e.g. P{framenr} gives P for framenr.
% 
% Note that "G" are ghost cells for which there's no match in a previous
% frame.
% 
% Hierarchy of calling:
% >> MW_makeSchnitzFileFromTracking 
%       >> DJK_data_treat
%       >> recalc_schnitz 
%           >> cellids

% Most comments by MW (added later)

%% Determine how many frame indices there are
highestFrameIdx = size(P,2);

%% The cells in the first frame can be assigned schnitz numbers starting from 1
nrCellsInFirstFrame = size(P{1},1);
cellid{1}=[1:nrCellsInFirstFrame]';

% Correspondingly, the schnitzes in the first frame do not have parents.
% (Not sure about DD & EE though).
%Make PP etc zero for numel(P{1})
PP(1:nrCellsInFirstFrame ) = 0;
DD(1:nrCellsInFirstFrame ) = 0;%ST
EE(1:nrCellsInFirstFrame ) = 0;%ST

% and record where our index should be now
nextid=nrCellsInFirstFrame;

%% Loop over remaining frames
% This is an iterative procedure, where if a cell did not divide, it'll
% be assigned the same schnitznr it had last frame. If it has divided, 
% we'll assign it a new schnitznr, and the next frame will inherit those
% numbers.
% Note that a dividing cell is encoded as follows:
% 0 d 0 p2
% 0 0 e p2
% (Where p2 is the celllabel in the n+1 frame.)
% And a non-dividing cell as 
% p1 0 0 p2
for previousFrameIdx=1:highestFrameIdx
        
    currentFrameIdx=previousFrameIdx+1;
    
    % Unused code 
    k=1:size(P{previousFrameIdx},1);
    l=P{previousFrameIdx}>0;
    
    % detect cells that don't divide, and let them inherit the schnitznr
    % from the last frame
    nonZeroIdxP = P{previousFrameIdx}>0;
    cellid{currentFrameIdx}( P{previousFrameIdx}(nonZeroIdxP) ) = ...
        cellid{previousFrameIdx}(nonZeroIdxP);
    
    %% find cells that divided last frame 
    
    % First those encoded by D array
    if any(D{previousFrameIdx}>0)
        
        % Find them
        nonZeroIdxD = find(D{previousFrameIdx}>0);
        
        % Assign those cells in the current frame new schnitznrs
        cellid{currentFrameIdx}( D{previousFrameIdx}( nonZeroIdxD ) ) = ...
            [nextid+1:nextid+numel(nonZeroIdxD) ];
        
        % Create PP, where
        % PP(schnitzindx) = the schnitzindex of the parent 
        PP( nextid+1 : nextid + sum(nonZeroIdxD) ) = ...
            cellid{previousFrameIdx}(nonZeroIdxD) ;
        
        % Update the next "free" index.
        nextid = nextid + numel(nonZeroIdxD);
        
    end
    
    % Then those encoded by E array (code identical to above)
    if any(E{previousFrameIdx}>0)
        % Find them
        nonZeroIdxE = find(E{previousFrameIdx}>0);
        
        % Assign those cells in the current frame new schnitznrs
        cellid{currentFrameIdx}( E{previousFrameIdx}( nonZeroIdxE ) ) = ...
            [nextid+1:nextid+numel(nonZeroIdxE) ];
        
        % Create PP, where
        % PP(schnitzindx) = the schnitzindex of the parent 
        PP( nextid+1 : nextid + numel(nonZeroIdxE) ) = ...
            cellid{previousFrameIdx}(nonZeroIdxE) ;
        
        % Update the next "free" index.
        nextid = nextid + numel(nonZeroIdxE);
    end
   
    % Then those encoded by G
    % These are special cases, since these are cells which do not have a
    % parent. Again, same code as above, except that PP(idx) will be 0,
    % because that's flag for orphan cell lines.
    if any(G{currentFrameIdx})
        warning('Ghost cells detected (i.e. incomplete lineages).');
        nonZeroIdxG = find(G{previousFrameIdx}>0); % Find
        cellid{currentFrameIdx}( nonZeroIdxG ) =[nextid+1:nextid+ sum(nonZeroIdxG)]; % Assign new nr
        PP(  [nextid+1 : nextid + numel(nonZeroIdxG) ] ) = 0; % 0 = has no parent
        nextid=nextid+ numel(nonZeroIdxG); % update nxt free
    end

    % Now create DD (the "reverse" of PP), where
    % DD(schnitzindx) = the schnitzindex of the daughter cell. 
    if any(D{previousFrameIdx}>0) 
        
        % Again, find them
        nonZeroIdxDprevious= find(D{previousFrameIdx}>0);
        
        % Create DD similar to PP
        DD(cellid{previousFrameIdx}(nonZeroIdxDprevious)) = ...
            cellid{currentFrameIdx}( D{previousFrameIdx}( nonZeroIdxDprevious ) );
    end
    
    % Now same for EE, which is equivalent to DD.
    if any(E{previousFrameIdx}>0)
        % Again, find them
        nonZeroIdxEprevious = find(E{previousFrameIdx}>0);
        % Create EE similar to PP
        EE(cellid{previousFrameIdx}( nonZeroIdxEprevious)) = ...
            cellid{currentFrameIdx}( E{previousFrameIdx}( nonZeroIdxEprevious ) );
    end
    
end


if( numel(DD) <nextid )
    DD(nextid)=0;
end

if( numel(EE) <nextid )
    EE(nextid)=0;
end

if( numel(PP) <nextid )
    PP(nextid)=0;
end















