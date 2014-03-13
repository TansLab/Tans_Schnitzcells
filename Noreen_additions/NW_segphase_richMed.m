% largely modified by Noreen Walker 2014-01
% original by Philippe Nghe 16/01/2012
%steps of the segmentation :
%   A. Find a mask
%   B. Find edges
%   C. Get clean cell mask and cut long cells
%   D. Remove fuzzy outline. Connect cut-off cell poles 
%   E. Create skeleton (seeds)
%   Z. Watershed from seeds
% INPUTS:
%   imageToSegment:    greyscale image
%   varargin:      parameters for image treatment indicated below
% OUTPUTS:
%   A_cropPhImage:  phase contrast image cropped (default is median filtered input image)
%   Z_segmentedImage:    segmented image (smaller size)
%   ROI_segmentation:   corner coordinates of the crop region in the original image

function [A_cropPhImage, Z_segmentedImage, ROI_segmentation] = NW_segphase_v1(imageToSegment,varargin)


%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%
%Reads parameter values or initializes with default values not too bad for
%E. Coli colonies on acryl gel pad with 100X phase contrast images.
q = inputParser;
q.addRequired('imageToSegment', @isnumeric);

%STEP A : finds a global mask and crop the image
q.addParamValue('rangeFiltSize',35,@isnumeric);       %dectection of interesting parts of the image
q.addParamValue('maskMargin',20,@isnumeric);          %additional margin of the mask : enlarge if cells missing on the edge

%STEP XYZ
q.addParamValue('LoG_Smoothing',2,@isnumeric);         %smoothing amplitude of the edge detection filter
q.addParamValue('minCellArea',250,@isnumeric);        %minimum cell area (smaller are erased)

%STEP XYZ
q.addParamValue('GaussianFilter',5,@isnumeric);        %smoothing of the original image to find local minima within cells 
q.addParamValue('minDepth',5,@isnumeric);             %minimum accepted depth for a local minimum

%STEP C : treatment of long cells
q.addParamValue('neckDepth',2,@isnumeric);    

%saving images
q.addParamValue('saveSteps',false,@islogical);  %indicate if you want to save intermediate images  
q.addParamValue('saveDir',pwd,@ischar);    %directory where intermediate treatment images are saved      

q.parse(imageToSegment, varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP O (O the letter, not 0 the zero !): suppress shot noise
O_PhImageFilt = medfilt2(imageToSegment);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP A : find a global mask and crop the image

A_maskImage = rangefilt(O_PhImageFilt,true(q.Results.rangeFiltSize));              %detect zones of sufficient intensity variations
A_maskImage = im2bw(A_maskImage,graythresh(A_maskImage));                  %threshold to black and white
A_maskImage = imclose(A_maskImage,strel('disk',q.Results.maskMargin));                     %enlarge mask         
labelMaskImage = logical(A_maskImage);                                     %only keep the biggest connected part of the mask
propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');
[~,idx] = max([propsMaskImage.Area]);
A_cropMaskImage = propsMaskImage(idx).Image;                               %cropped mask
bb = propsMaskImage(idx).BoundingBox;
ROI_segmentation = [ceil(bb([2 1])) ceil(bb([2 1]))+bb([4 3])-[1 1]];   %ymin xmin ymax xmax
A_cropPhImage = O_PhImageFilt(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4)); %cropped ph image image

if q.Results.saveSteps
    savePNGofImage(A_maskImage,'A_maskImage',q.Results.saveDir);
    %saveas(A_maskImage,'A_maskImage','png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP B : find edges 

B_negPh = imcomplement(A_cropPhImage);
se = strel('disk',1);
B_negPhErode = imerode(B_negPh, se);                                       %Morphological reconstruction ...
B_negPh = imreconstruct(B_negPhErode,B_negPh);                             %... in the mask of the original image
B_negPh = imdilate(B_negPh, se);                                          %Allows a better edge determination at contacts.
B_edgeImage1 = edge(B_negPh,'log',0,q.Results.LoG_Smoothing);                      %edge detection by smoothing (laplacian of gaussian ensures closed edges)

%suppress noisy surroundings
B_edgeImage2 = B_edgeImage1 & A_cropMaskImage;
B_fillEdgeImage2 = imfill(B_edgeImage2,'holes');
B_fillEdgeImage2 = bwareaopen(B_fillEdgeImage2,q.Results.minCellArea,4);           %suppress small stuff 
B_edgeImage2 = bwareaopen(B_edgeImage2,30);           %obsolete? blubb                     %remove small loops related to intracell variations


if q.Results.saveSteps
savePNGofImage(B_edgeImage1,'B_edgeImage1',q.Results.saveDir);
savePNGofImage(B_edgeImage2,'B_edgeImage2',q.Results.saveDir);
savePNGofImage(B_fillEdgeImage2,'B_fillEdgeImage2',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP C : get clean cell mask (dirt removed) and cut long cells

%makes a clean filled cells mask (removes segmented dirt without clear
%intensity-min)
C_smoothPh = imfilter(A_cropPhImage,fspecial('gaussian',q.Results.GaussianFilter,q.Results.GaussianFilter));
C_localMinPh = imextendedmin(C_smoothPh,q.Results.minDepth) & B_fillEdgeImage2;       %local minima in the mask      
C_cellMask = imfill(B_edgeImage2,find(C_localMinPh));                       
C_cellMask = bwmorph(C_cellMask,'open');                                    %clean cells mask

%shrinking steps to cut some cells
C_cellMaskSmall = C_cellMask & ~B_edgeImage2;                                      %thins cells by removing the edge

%cut long cells
continueToCut = true;
icut=bwmorph(C_cellMaskSmall,'thin',inf);

brchptsicut = PN_FindBranchPoints(icut);
icut(logical(brchptsicut)) = 0;
%some cleaning
icut = bwmorph(icut,'spur',3);
icut = bwareaopen(icut,10,8);

while continueToCut
    [cellsToRemove cutPoints] = NW_CutLongCells_richMed(icut,C_cellMaskSmall,q.Results.neckDepth);
    if max(max(cutPoints))==0
        continueToCut = false;
    else
        cutPoints = bwmorph(cutPoints,'dilate',2);
        C_cellMaskSmall(cutPoints)=0; %it could happen that not enough pixels are removed! but larger dilation -> stranger cell-pole shape
    %    C_seeds2(cutPoints) = false; %cuts the long cells on the seeds image
        icut(cutPoints) = false;
        icut= icut & ~cellsToRemove;
    end
end

if q.Results.saveSteps
    savePNGofImage(C_cellMask & ~B_edgeImage2,'C_cellMaskSmall_preCutLong',q.Results.saveDir);
    savePNGofImage(C_cellMaskSmall,'C_cellMaskSmall',q.Results.saveDir);
    savePNGofImage(C_cellMaskSmall & ~C_localMinPh,'C_cellMaskSmall and minima',q.Results.saveDir);
end
disp('finished *C* clean cellMask & cut Long Cells.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP D : Remove fuzzy outline & gulfs. Connect cut-off cell poles to
% larger cells. specific for rich medium

% ------ 
% connect suspiciously small cell areas (cutoff-poles...) to neighbour
% cells
% label all areas
D_cellMask_conn=bwlabel(C_cellMaskSmall);

areacutoff=400; % adjust
allAreas=regionprops(D_cellMask_conn,'Area');
idxSmallArea=find([allAreas.Area]<areacutoff);
for i=1:length(idxSmallArea);
    myareaidx=idxSmallArea(i);
    % check if area is still suspiciously small or already connected
    currentarea=regionprops(D_cellMask_conn==myareaidx,'Area');
    if isempty(currentarea) | currentarea.Area>areacutoff; continue; end
    D_cellMask_conn=NW_ConnectDividedCellParts(D_cellMask_conn,B_negPh,myareaidx);
end

% -----
% imclose each area
mystruct_el=strel('disk',10); %adjust
D_cellMask_closed=NW_imclose_eachArea(D_cellMask_conn,mystruct_el);
% reduce image again to logical: areas=1, background=0
D_cellMask_closed=(D_cellMask_closed>0);

E_seeds = bwmorph(D_cellMask_closed,'thin',inf);     %for later                     %skeletinization

% reconstruct Cell_Mask by dilation. compensate for previous edge removal
D_cellMask_closed=imdilate(D_cellMask_closed,strel('disk',1));

if q.Results.saveSteps
    savePNGofImage(D_cellMask_conn,'D_cellMask_conn',q.Results.saveDir);
    savePNGofImage(D_cellMask_closed,'D_cellMask_closed',q.Results.saveDir);
end
disp('Finished *D* imclose and connect small cells.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP E : create skeleton and suppress branch points of the skeleton
%E_seeds = bwmorph(D_cellMask_closed,'thin',inf);                                    %skeletinization
brchpts = PN_FindBranchPoints(E_seeds);
E_seeds(logical(brchpts)) = 0;
%some cleaning
E_seeds = bwmorph(E_seeds,'spur',3);
E_seeds = bwareaopen(E_seeds,10,8);

if q.Results.saveSteps
    savePNGofImage(E_seeds,'E_seeds',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP Z : final segmentation by watershedding

%prepare mask for watershedding
Z_maskToFill = bwmorph(D_cellMask_closed,'dilate');
Z_background = ~Z_maskToFill;                                               %background

%prepare final seeds
Z_seeds = bwmorph(E_seeds,'spur',3);
Z_seeds = bwareaopen(Z_seeds,4,8);

%prepare seeded landscape for watershedding
Z_d1 = -bwdist(Z_background);                                               %distance transform  
Z_d1 = imimposemin(Z_d1, Z_background | Z_seeds);                           %impose background and seeds

%watershedding
Z_segmentedImage = watershed(Z_d1);                     
Z_segmentedImage(Z_background) = 0;                                         % (small) areas which did not contain seed are now also included
Z_segmentedImage = bwareaopen(Z_segmentedImage,10);
Z_segmentedImage = bwlabel(Z_segmentedImage);                               %segmentation is finished here
Z_segmentedImage = imdilate(Z_segmentedImage, strel('diamond',1));

if q.Results.saveSteps
savePNGofImage(Z_maskToFill & ~Z_seeds,'Z_Mask and seeds',q.Results.saveDir);
savePNGofImage(Z_segmentedImage,'Z_segmentedImage',q.Results.saveDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION END   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
if max(max(Z_segmentedImage))==0                                           %let user know when no cells were found
  disp([' * !! WATCH OUT !! no cells found on this frame...']);
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function savePNGofImage(image, name, saveDirectory)
    filename = [saveDirectory,name,'.png'];
    if isempty(image)
        image = [ 0 ]; % will get error if image is empty
    end
    DJK_writeSegImage(image, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
