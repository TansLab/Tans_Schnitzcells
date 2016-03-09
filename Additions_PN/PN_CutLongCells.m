%by Philippe Nghe 16/01/2012

function [nonCutCells cutImage] = PN_CutLongCells(skelim,refim,neckDepth)
% Parameters that can be set in this function: 
%
% LONGFACTOR:       determines when to inspect cells for division; based on
%                   ratio to mean length.
% MINCELLS :        number of cells below which absolute cutoff for
%                   identifying suspiciously long (possibly dividing) is
%                   used.
% SUSPICIOUSSIZE :  the absolute size mentioned above.

% CONFIG VARS FOR THIS FUNCTION--------------------------------------------
LONGFACTOR = 1.5;
MINCELLS = 8;   % MW added to handle cuts of very small colonies, below 
                % this number of cells, absolute length will be used
SUSPICIOUSSIZE = 30; % If #cells < mincells, use this criterion to 
                     % detect long cells
% -------------------------------------------------------------------------



sref = size(refim);
cutImage = zeros(sref);
nonCutCells = zeros(sref);

dimage = bwdist(~refim);%distance transform
cc=bwconncomp(skelim);
stats = regionprops(cc,'Area','BoundingBox','Image');
characSize = median([stats.Area]);

% detect cells that are long relative to the mean of all cells, 
% these cells will be analyzed for cell division
idx = find([stats.Area] > characSize*LONGFACTOR);  %suspicious cells
                                            % parameter can be changed for 
                                            % different conditions.

% MW addition 2014/06/22
% If the colony has few cells, deviation from average length doesn't tell 
% you much, and cells might not be cut. Following code in that situation
% adds cells with a large absolute size (defined by SUSPICIOUSSIZE) to the
% index of cells to be inspected.
if length(stats)<MINCELLS
   idx = unique([idx find([stats.Area] > SUSPICIOUSSIZE)]);
   
   if isempty(idx) % this is to deal with the fact that unique([]) returns 
                   % an empty 0-by-1 matrix, which will start loop below..
                   % MW 2014/11
       idx = [];
   end
end

for ii = idx  %study the case of long cells individually

    %extracts the subimage containing the long cell
    s = stats(ii);
    xb = ceil(s.BoundingBox(1));    yb = ceil(s.BoundingBox(2));
    lx = s.BoundingBox(3);          ly = s.BoundingBox(4);
    subImage = imcrop(dimage,[xb yb lx-1 ly-1]);
    localSkel = s.Image;    %local skeleton
    subImage(~localSkel) = inf; %local restriction of distance transform to the skeleton 
    
    %cuts under the condition that the necking has a certain depth
    [m xm ym] = MinCoordinates2(subImage);    %potentiel division point
    localSkel(ym,xm)=false; %cut...
    localCc=bwconncomp(localSkel); %...and examin the 2 sides of the cells
    cutHappened = false;
    if localCc.NumObjects == 2 %avoid cases in which the cut point is at a cell end
        av_left = mean(subImage(localCc.PixelIdxList{1}));    %average thickness on one side
        av_right = mean(subImage(localCc.PixelIdxList{2}));    %average thickness on other side
        if (av_left-m > neckDepth) && (av_right-m > neckDepth) %cusp of sufficient depth
            cutImage(ym+yb-1,xm+xb-1)=1;
            cutHappened = true;
        end
    end
    
    %cells which have not been cut are stored for erasion to not be reexamined
    if ~cutHappened
       nonCutCells(cc.PixelIdxList{ii}) = 1; 
    end
end

clear sref dimage cc stats characSize idx s xb yb lx ly m xm ym localSkel localCc cutHappened av_left av_right
end