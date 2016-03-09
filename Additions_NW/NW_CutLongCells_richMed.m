%modified by Noreen Walker 2013-01 to account for segmentation issues in
%rich medium (substructures)
%Original by Philippe Nghe 16/01/2012

function [nonCutCells cutImage] = NW_cutLongCells_richMed(skelim,refim,neckDepth)

% **** ADJUST ****
suspiciousSizeFactor=2;   % use=2 (contrary to minimal Medium=1.5)
remEdgePx=10;%20;             % remove xx end pixels of skeletons where cut must not happen
% ****************

sref = size(refim);
cutImage = zeros(sref);
nonCutCells = zeros(sref);

dimage = bwdist(~refim);%distance transform
cc=bwconncomp(skelim);
stats = regionprops(cc,'Area','BoundingBox','Image');
characSize = median([stats.Area]);
idx = find([stats.Area] > characSize*suspiciousSizeFactor); %suspicious cells 

for ii = idx  %study the case of long cells individually
    
    %extracts the subimage containing the long cell
    s = stats(ii);
    xb = ceil(s.BoundingBox(1));    yb = ceil(s.BoundingBox(2));
    lx = s.BoundingBox(3);          ly = s.BoundingBox(4);
    subImage = imcrop(dimage,[xb yb lx-1 ly-1]);
    localSkel = s.Image;    %local skeleton
    localSkel=bwmorph(localSkel,'shrink',remEdgePx); % remove end pixels (rich Med specific)
    subImage(~localSkel) = inf; %local restriction of distance transform to the skeleton 
    
    %cuts under the condition that the necking has a certain depth

    % Tests several potential neck positions (specific for rich medium ->
    % branched skeletons)
    if     ii==158
        disp(' ')
    end
    % [m xm ym] = MinCoordinates2(subImage);    %potentiel division point
    [mvec xmvec ymvec] = MinCoordinates2_N(subImage,3);    %potentiel division point, 3 closest dist's
    
    % loop over potential division points
    for potdivrun=1:length(xmvec)
        m=mvec(potdivrun);
        xm=xmvec(potdivrun);
        ym=ymvec(potdivrun);
        localSkelRun=localSkel;
        localSkelRun(ym,xm)=false; %cut...
        localCc=bwconncomp(localSkelRun); %...and examin the 2 sides of the cells
        cutHappened = false;
        if localCc.NumObjects == 2 %avoid cases in which the cut point is at a cell end
            av_left = mean(subImage(localCc.PixelIdxList{1}));    %average thickness on one side
            av_right = mean(subImage(localCc.PixelIdxList{2}));    %average thickness on other side
            if (av_left-m > neckDepth) && (av_right-m > neckDepth) %cusp of sufficient depth
                cutImage(ym+yb-1,xm+xb-1)=1;
                cutHappened = true;
                break %stop when cell could be cut. maybe allow for more cuts in future
            end

        end
    end
    
    %cells which have not been cut are stored for erasion to not be reexamined
    if ~cutHappened
       nonCutCells(cc.PixelIdxList{ii}) = 1; 
    end
end

clear sref dimage cc stats characSize idx s xb yb lx ly m xm ym localSkel localCc cutHappened av_left av_right
end