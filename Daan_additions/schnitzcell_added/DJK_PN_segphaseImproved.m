function [phsub, L9C, rect, phsegsub] = DJK_segphaseImproved(ph3, p)
% function [phsub, L9C, rect, phsegsub] = DJK_segphaseImproved(ph3, p)
% 
% copied and rewritten from DJK_segphase, which was copied from segphase
%
% INPUTS:
%   ph3:    phase contrast images (multiple Z slices)
%   p:      schnitzcells settings
%
% OUTPUTS:
%   phsub:  phase contrast images (smaller size, scaled constrast)
%   L9C:    segmented image (smaller size)
%   rect:   transformation required to reconstruct full size image from smaller ones
%   phseg:  phase contrast image used for edges (default is mean of ph3)
%
% SAVES IMAGES:
%   L:      segmented image (full size)
%
%
% SEGMENTATION STEPS:
%   1. finds edge
%   2. finds mask
%   3. fills cells using imcomplement
%   4. fills edge using cells found from incomplement
%   5. breaks up cell clumps (cutcurvbpts)
%   6. cuts up individual cells (cutcurv)
%   7. cuts up any long cells remaining (breakcell)
%   8. cuts up any kinky cells remaining (breakcell, dekinker)
%   9. removes any small `cells' left
%  10. relabels image

%--------------------------------------------------------------------------
% SETTINGS
%--------------------------------------------------------------------------
p.numphaseslices; % number of slices
image_size = size(ph3(:,:,1)); % size of original image
image_center = image_size/2; % center of image

% STEP 0A: FILTER AND EDGE PHASE
  p.edgeSlices; % slices that are used for edge detection
  p.useMedfilt2forEdge; % whether medfilt2 is applied on images for edge (default 0)
  p.edge_lapofgauss_sigma; % important for edge detection. Higher gives more smoother edges

% STEP 0B: FIND MASK_1 (smallest region containing cells)
  boxsize= 40; % size of box that is run over the edge image for numEdgePixel counting
  replaceboxsize = 5; % Step size for box (boxsize is shrunk to replaceboxsize on output image)
  p.minNumEdgePixels; % mask will be where less edge pixels in box than p.minNumEdgePixels
  % imopen with strel('disk',25)

% DELETED: STEP 0C: REDUCE MASK_1 (only keep large clumps)

% STEP 0D: EXTRACT SMALLER SUBREGIONS
  extraPixels = 5; % extra number of pixels on either side of region with segmented cells

% STEP 0E: MAKE AN EVEN SMALLER MASK: MASK_2
  % imclose with strel('disk',9)
  % imdilate with strel('disk',5)
  % imdilate with strel('disk',15)

% STEP 1: FILL CELLS 
  % imdilate with strel('disk',10)
  p.L1Fill_only_use_mask_2_bigger; % option to only use mask_2_bigger, and not mask_1 (default: 0)

% STEP 2: USING BOTTOM HAT, GET SEEDS TO REMOVE REGION BETWEEN CELLS
  p.fillingEdgeSlices; % slices that are used for filling the edges with bottom hat
  p.botHatSize; % imbothat with strel('disk', p.botHatSize)

% STEP 3: CLEAN UP UNSIGHTLY CELLS
  p.minCellLengthConservative; % remove cells where MajorAxisLength < p.minCellLengthConservative

% STEP 4: FILL IN EDGE (using seeds from L3)
  strelsize = p.botHatSize - 1; % reduce result from botHat with this many pixels
  % imdilate with strel('disk',3) % remove stuff 3 pixels outside mask_2

% STEP 5A: CUTTING BY REMOVING SINGLE PIXELS THAT CONNECT DISPARATE CELLS

% STEP 5B: CUTTING AT BRANCHPOINTS USING cutcurvbpts
  p.maxThreshCut; % 
  p.maxCellWidth; % 

% STEP 5C: REMOVE SMALL CELLS
  % removes cells where Area < 5

% STEP 6A: CUTTING CELLS AT CONCAVE POINTS USING cutcurv
  % cuts cells where Solidity < 0.85
  p.maxThreshCut2; % smaller maxThreshCut2 the more points are cut 
  p.maxCellWidthConservative; % two points must be closer than this for cut to be accepted
  p.minCellLength; % cut ignored if it creates `cell' smaller than this

% STEP 6B: CUTTING CELLS AT CONCAVE POINTS USING cutcurv
  % cuts cells where Solidity < 0.85
  p.maxThreshCut2; % smaller maxThreshCut2 the more points are cut 
  p.maxCellWidth; % two points must be closer than this for cut to be accepted
  p.minCellLength; % cut ignored if it creates `cell' smaller than this

% STEP 7: CUT LONG CELLS BY LOOKING AT PHASECONTRAST USING breakcell
  % cut cells where MajorAxisLength > 70 
  p.edgeSlices; % slices that are used to look for changes in phase along thin
  p.maxThresh; %
  p.minThresh; %
  p.minCellLength; % cut ignored if it creates `cell' smaller than this

% DELETED: STEP 8: DEKINKING CELLS WITH LOW SOLIDITY

% STEP 9A: REMOVE SMALL CELLS
  p.minCellLengthConservative; % remove cells where MajorAxisLength < p.minCellLengthConservative
  
% STEP 9B: INCREASE CELL WITH 1 PIXEL SO EDGE PIXELS ARE INCLUDED
  % imdilate with strel('disk',1)

% STEP 9C: RESIZE OUTPUT IMAGES TO EXACTLY FIT SEGMENTED CELLS
  extraPixels; % uses extraPixels for extra pixels on sides
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 0A: FILTER AND EDGE PHASE
%--------------------------------------------------------------------------
% Original Process here:
%   * each slice in ph3 was medfilt2 with [3 3]
%   * each slice was rescaled, such that 25th highest value is 10000 and the
%     25th lowest value is 0
%   * edge was performed on p.segmentationPhaseSlice only
%
% New Process:
%   * if p.useMedfilt2forEdge is 1, medfilt2 is applied on slices in ph3
%   * phseg is mean of slices in ph3 indicated in p.edgeSlices 
%   * contrast of phseg and slices in ph3 are scaled
%   * edge is determined from phseg and p.edge_lapofgauss_sigma

% median filter 
if p.useMedfilt2forEdge
  for i = 1:size(ph3,3) % loop over slices
    ph3(:,:,i)= medfilt2(ph3(:,:,i),[3 3]);
  end
end

% make phseg
phseg = mean( ph3(:, :, p.edgeSlices), 3);

% rescale phase images
for i = 1:size(ph3,3) % loop over slices
  ph3(:,:,i) = scaleImage( ph3(:,:,i) );
end
phseg = scaleImage( phseg );

% edge detection
me = edge( phseg, 'log', 0, p.edge_lapofgauss_sigma);

savePNGofImage(me,'L0A_me',p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 0B: FIND MASK_1 (smallest region containing cells)
%--------------------------------------------------------------------------
% Process:
%   * Run box of size boxsize over the edge image in steps of replaceboxsize
%   * count the number of white pixels in the box, cells have less than minNumEdgePixels pixels
%   * tidy up and slightly enlarge mask (imdilate & imopen)

% mask_1 will be binary image where 0 means area with no cells, 1 area with cells
% used to be imout
mask_1 = zeros(image_size);

% loop box over image
for x = 1:replaceboxsize:image_size(1) - boxsize
  for y = 1:replaceboxsize:image_size(2) - boxsize
    xmin = x + boxsize/2;
    xmax = x + boxsize/2 + replaceboxsize - 1;
    ymin = y + boxsize/2;
    ymax = y + boxsize/2 + replaceboxsize - 1;
    subregion = me(x:x + boxsize, y:y + boxsize);
    checkiscell = sum(subregion(:)) < p.minNumEdgePixels;
    mask_1(xmin:xmax, ymin:ymax) = ones(replaceboxsize)*checkiscell;
  end
end

% dilate output image to make sure ends of cells included
mask_1 = uint8(imdilate(mask_1,strel('disk',replaceboxsize)));

% tidy up imout so that only one highlighted region remains
mask_1_temp = imopen(mask_1,strel('disk',25));
if max(mask_1_temp(:))
    mask_1 = mask_1_temp;
end

savePNGofImage(mask_1,'L0B_mask_1',p); 
clear mask_1_temp;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 0D: EXTRACT SMALLER SUBREGIONS
%--------------------------------------------------------------------------
% Process:
%   * based on mask_1, extract smaller subregions of used images
%   * extraPixels number of extra pixels that are added on each side

% Get bounding box (rect) of mask_1 with extraPixels
[fx,fy] = find(mask_1);
xmin = max(min(fx) - extraPixels, 1);
xmax = min(max(fx) + extraPixels, size(mask_1,1));
ymin = max(min(fy) - extraPixels, 1);
ymax = min(max(fy) + extraPixels, size(mask_1,2));
rect = [xmin ymin xmax ymax];

% extract smaller subregions
phsub = ph3(xmin:xmax, ymin:ymax, :);
phsegsub = phseg(xmin:xmax, ymin:ymax);
mask_1sub = mask_1(xmin:xmax, ymin:ymax); % DJK: was imout2
mesub = me(xmin:xmax, ymin:ymax);

% savePNGofImage(mask_1sub,'L0D_mask_1sub',p);
% savePNGofImage(mesub,'L0D_mesub',p);
clear fx fy xmin xmax ymin ymax;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 0E: MAKE AN EVEN SMALLER MASK: MASK_2
%--------------------------------------------------------------------------
% Process:
%   * get another mask by using a threshold
%   * use imclose, imdilate, imfill to improve
%   * make a bigger version by imdilate

% use a threshold to find where the cells are
bw = im2bw( phsegsub, joeslowgraythresh(phsegsub) );
wb = imcomplement(bw);

% smooth edges
wbc = imclose(wb, strel('disk',9));

% JCR: add back thresholded cells that were lost by the close, for edge issues
wbc = wbc | wb;

% fill and dilate
mask_2 = imdilate(imfill(wbc,'holes'), strel('disk',5));
mask_2 = imfill(mask_2,'holes');

% JCR trying larger dilation to handle fuzzy cells at edge of colony
% JCR increasing to avoid losing entire colonies (when colony gets bright)
mask_2_bigger = imdilate(mask_2, strel('disk',10)); % DJK 081228 used to be 15

% JCR adding extra hole filling after second dilation
mask_2_bigger = imfill(mask_2_bigger,'holes');

savePNGofImage(mask_2,'L0E_mask_2',p);
savePNGofImage(mask_2_bigger,'L0E_mask_2_bigger',p);
clear bw wb wbc;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 1: FILL CELLS 
%--------------------------------------------------------------------------
% Process:
%   * get edge within both mask_1 & mask_2_bigger
%   * fill all regions not open to background

% get masked part of mesub
if p.L1Fill_only_use_mask_2_bigger
  mesub_mask = mask_2_bigger & mesub;
else
  mesub_mask = mask_2_bigger & mesub & imdilate(mask_1sub,strel('disk',10));
end

% Give each region a number
L1 = bwlabel(imcomplement(mesub_mask), 4);

% Give background value 0
r = regionprops(L1,'Area');
backgroundObject = find([r.Area]==max([r.Area]));
L1temp = L1;
L1temp(L1==backgroundObject) = 0;
L1temp(L1==0) = backgroundObject;
L1 = bwlabel(L1temp,4);

savePNGofImage(mesub_mask,['L1_mesub_mask'],p);
savePNGofImage(L1,'L1',p);
clear r backgroundObject L1temp;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 2: GET SEEDS TO REMOVE REGION BETWEEN CELLS
%--------------------------------------------------------------------------

L2 = L1;

% uses the image that is the mean of p.fillingEdgeSlices
phsegFill = mean( phsub(:, :, p.fillingEdgeSlices), 3);

%finds local minima in the mask
phsegFillgaussian = imfilter(phsegFill,fspecial('gaussian',5,5));                
localMin = imhmin(phsegFillgaussian,5) ; 
L2 = localMin & imdilate(mask_2_bigger, strel('disk',p.L2_increase_mask_2_bigger));
[fx,fy] = imregionalmin(L2);   %seeds

savePNGofImage(L2,'L2',p);
clear phsegFill phsegFillgaussian localMin;
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% STEP 4: FILL IN EDGE (using seeds from L3)
%--------------------------------------------------------------------------

  L4 = imfill(mesub,[fx,fy]);

% Label cells
L4 = bwlabel(L4);
L4 = renumberimage(L4); 

savePNGofImage(L4,'L4',p);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 5A: CUTTING BY REMOVING SINGLE PIXELS THAT CONNECT DISPARATE CELLS
%--------------------------------------------------------------------------
% Process:
%   * removes single pixels that connect otherwise disparate regions

% rmsinglepointconnections removes single pixels that connect otherwise
% disparate regions of image, after infinite bwmorph thinning
L5A = rmsinglepointconnections(L4);

savePNGofImage(L5A,'L5A',p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 5B: CUTTING AT BRANCHPOINTS USING cutcurvbpts
%--------------------------------------------------------------------------
% Process:
%   * morphologically thin the cell-segment to a single-pixel line
%   * find branchpoints in cell, i.e. points where several cells are grouped together in a clump
%   * cut branchpoints

L5B = L5A;

for i = 1:max2(L5B)
  % get only the cell that is gonna be cut
  Lcell = (L5B == i);

  % try to cut cell
  cutcell = DJK_cutcurvbpts(Lcell, p.maxThreshCut, p.maxCellWidth); % ,1 for figs
 
  % remove original cell
  L5B(Lcell) = 0;

  % place cutcell (can be 1 or more cells)
  cellnos = unique(cutcell);  
  label = max2(L5B);
  for j = 2:length(cellnos)
    L5B(find(cutcell == cellnos(j)))= label+j;
  end
end
L5B = renumberimage(L5B);

savePNGofImage(L5B,'L5B',p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 5C: REMOVE SMALL CELLS
%--------------------------------------------------------------------------
% Process:
%   * remove cells with area < 5

L5C = L5B;

% Removes cells that are just small pixels
r = regionprops(L5C, 'area');
fpts = find([r.Area]<5);
for i = 1:length(fpts),
  L5C(L5C == fpts(i))= 0;
end;
L5C = renumberimage(L5C);

savePNGofImage(L5C,'L5C',p);
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% STEP 6A: CUTTING CELLS AT CONCAVE POINTS USING cutcurv
%--------------------------------------------------------------------------
% Process:
%   * only performed if p.waistCut = 1 (default is 0)
%   * find locations where solidity < 0.85
%   * try to cut there using cutcurv, width must be > maxCellWidthConservative & result in cells > minCellLength
%   
% DJK: Added before 6B, so first cuts are applied at narrow locations 
%      (uses maxCellWidthConservative in stead of maxCellWidth)
%
% Further info: cuts cells at the narrow waist: makes sure septation events
% are properly identified and that cells are identified seperately. Cuts
% cells at points where both sides of the cell are sufficiently concave.

L6A = L5C;

if p.waistCut & p.maxCellWidthConservative < p.maxCellWidth
  r = regionprops(L5C, 'solidity');
  fkinks = find(([r.Solidity] < 0.85));
  for i = 1:length(fkinks)
    % get only the cell that is gonna be cut
    Lcell = (L5C == fkinks(i));
  
    % try to cut cell
    cutcell = DJK_cutcurv(Lcell, p.maxThreshCut2, p.maxCellWidthConservative, p.minCellLength);

    % remove original cell
    L6A(Lcell) = 0;

    % place cutcell (can be 1 or more cells)
    cellnos = unique(cutcell);  
    label = max2(L6A);
    for j = 2:length(cellnos)
      L6A(find(cutcell == cellnos(j)))= label+j;
    end
  end
  L6A = renumberimage(L6A);
  
  savePNGofImage(L6A,'L6A',p);
end
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% STEP 6B: CUTTING CELLS AT CONCAVE POINTS USING cutcurv
%--------------------------------------------------------------------------
% Process:
%   * same as 6A, but cuts are allowed at locations where cell is more wide 

L6B = L6A;

if p.waistCut
  r = regionprops(L6A, 'solidity');
  fkinks = find(([r.Solidity] < 0.85));
  for i = 1:length(fkinks)
    % get only the cell that is gonna be cut
    Lcell = (L6A == fkinks(i));

    % try to cut cell
    cutcell = DJK_cutcurv(Lcell, p.maxThreshCut2, p.maxCellWidth, p.minCellLength);

    % remove original cell
    L6B(Lcell) = 0;

    % place cutcell (can be 1 or more cells)
    label = max2(L6B);
    cellnos = unique(cutcell);  
    for j = 2:length(cellnos)
      L6B(find(cutcell == cellnos(j)))= label+j;
    end
  end
  L6B = renumberimage(L6B);
  
  savePNGofImage(L6B,'L6B',p);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 7: CUT LONG CELLS BY LOOKING AT PHASECONTRAST USING breakcell
%--------------------------------------------------------------------------
% Process:
%   * find cells whose MajorAxisLength > 
%   * try to cut where there is a change in the phase value using phseg created from p.edgeSlices
%   * must result in cells > minCellLength

L7 = L6B;

r = regionprops( L6B, 'majoraxislength' );
fbiggies = find(([r.MajorAxisLength]>70)); 
for i = 1:length(fbiggies)
  % get only the cell that is gonna be broken
  Lcell = +(L6B == fbiggies(i)); % + converts logical to double
  Lcell(Lcell == 1)= fbiggies(i);
  % savePNGofImage(Lcell,['L7_' num2str(i)],p);
  
  % try to break cell, DJK 081229 used to use phsub(:,:,p.imNumber1), now phsegsub
  cutcell = breakcell(Lcell, phsegsub, p.maxThresh, p.minThresh, p.minCellLength);
  
  % remove original cell
  L7(L6B == fbiggies(i))= 0;
  
  % place cutcell (can be 1 or more cells)
  cellnos = unique(cutcell);  
  label = max2(L7);
  for j = 2:length(cellnos)
    L7(find(cutcell == cellnos(j)))= label+j;
  end
end
L7 = renumberimage(L7);

savePNGofImage(L7,'L7',p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 8: DELETED
%--------------------------------------------------------------------------
% DJK: this step works badly with bended cells, as these have low solidity.
%      so I deleted it

L8B = L7;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 9A: REMOVE SMALL CELLS
%--------------------------------------------------------------------------
% Process:
%   * remove cells with MajorAxisLength < p.minCellLengthConservative

L9A = L8B;

if p.removeSmallCells
  r = regionprops( L9A, 'majoraxislength' );
  flittles = find( [r.MajorAxisLength] < p.minCellLengthConservative );
  for i = 1:length(flittles)
    L9A(L9A == flittles(i))= 0; % delete cell
  end
  L9A = renumberimage(L9A);

  savePNGofImage(L9A,'L9A',p);
  clear r flittles;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 9B: INCREASE CELL WITH 1 PIXEL SO EDGE PIXELS ARE INCLUDED
%--------------------------------------------------------------------------
% Process:
%   * imdilate with 1 pixel

% DJK: old was original, Sander increased to 2 pixels, cause 'looks nicer'
L9B = imdilate( L9A, strel('diamond',1) );
L9B_old = carefuldilate(+L9A, strel('diamond',1), 1);
L9B_sander = carefuldilate(+L9A, strel('diamond',1), 2);

savePNGofImage(L9B,'L9B',p);
savePNGofImage(L9B_old,'L9B_old',p);
savePNGofImage(L9B_sander,'L9B_sander',p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% STEP 9C: RESIZE OUTPUT IMAGES TO EXACTLY FIT SEGMENTED CELLS
%--------------------------------------------------------------------------
% Process:
%   * get correct output for phsub, L9C, rect, phsegsub
%   * images are cropped to exactly fit segmented cells

L9B_fullsize = zeros( image_size );
L9B_fullsize( rect(1):rect(3), rect(2):rect(4) ) = L9B;
[fx,fy] = find(L9B_fullsize);
xmin = max(min(fx) - extraPixels, 1);
xmax = min(max(fx) + extraPixels, size(L9B_fullsize,1));
ymin = max(min(fy) - extraPixels, 1);
ymax = min(max(fy) + extraPixels, size(L9B_fullsize,2));

phsub = ph3( xmin:xmax, ymin:ymax, :);
L9C = L9B_fullsize(xmin:xmax, ymin:ymax);
rect = [xmin ymin xmax ymax];
phsegsub = phseg( xmin:xmax, ymin:ymax);

% saveTIFFofPhaseImage(phsegsub, 'phsegsub', p); % DJK 090530 turned off
savePNGofImage(L9C,'L9C',p);
%--------------------------------------------------------------------------

% AT THE END, just let user know when no cells were found
if max2(L9C)==0
  disp([' * !! WATCH OUT !! no cells found on this frame...']);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to scale image, such that 25th highest value is 10000 and the
% 25th lowest value is 0
function imout = scaleImage(imin)
x = double(imin);
s = sort(x(:));
small = s(25);
big = s(end-25);
rescaled = (x - small)/(big - small);
rescaled(rescaled<0) = 0;
imout = uint16(10000*rescaled);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function savePNGofImage(image, name, p)
  if isfield(p,'DJK_saveDir')
    filename = [p.segmentationDir,p.DJK_saveDir,name,'.png'];
    if isempty(image)
      image = [ 0 ]; % will get error if image is empty
    end
    DJK_writeSegImage(image, filename);
    disp(['          * written : ' name '.png']);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function saveTIFFofPhaseImage(image, name, p)
  if isfield(p,'DJK_saveDir')
    filename = [p.segmentationDir,p.DJK_saveDir,name];
    if isempty(image)
      image = [ 0 ]; % will get error if image is empty
    end
    DJK_saveTIFF(uint16(image), [p.segmentationDir,p.DJK_saveDir,name]);
    disp(['          * written : ' name '.tif']);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%