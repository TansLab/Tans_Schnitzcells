



% FUNCTION IS DEPRECATED
%
% Use option p1.useFullImage = 1 with PN_segmoviephase_3colors instead.
% This calls the PN_segphase.m function with the right options.
% 
% This function will be removed.
% - MW 2015/04






% NW (2013): No restriction to largest connected area -> cells that are far
% apart are all detected. Segmentation is not restricted to central area.
% Otherwise identical to Std Segmentation
%
%by Philippe Nghe 16/01/2012
%steps of the segmentation :
%   A. Find a mask
%   B. Find edges
%   C. Find skeletenized seeds
%   D. Cut branch points
%   E. Cut long cells
%   Z. Watershed from seeds
% INPUTS:
%   imageToSegment:    greyscale image
%   varargin:      parameters for image treatment indicated below
% OUTPUTS:
%   A_cropPhImage:  phase contrast image cropped (default is median filtered input image)
%   Z_segmentedImage:    segmented image (smaller size)
%   ROI_segmentation:   corner coordinates of the crop region in the original image

function [A_cropPhImage, Z_segmentedImage, ROI_segmentation] = PN_segphase(imageToSegment,varargin)

disp('WARNING: This function is deprecated, use p1.useFullImage = 1 with PN_segmoviephase_3colors instead.');

%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%
%Reads parameter values or initializes with default values not too bad for
%E. Coli colonies on acryl gel pad with 100X phase contrast images.
q = inputParser;
q.addRequired('imageToSegment', @isnumeric);

%STEP A : finds a global mask and crop the image
q.addParamValue('rangeFiltSize',35,@isnumeric);       %dectection of interesting parts of the image
q.addParamValue('maskMargin',20,@isnumeric);          %additional margin of the mask : enlarge if cells missing on the edge

%STEP B : find edges
q.addParamValue('LoG_Smoothing',2,@isnumeric);         %smoothing amplitude of the edge detection filter
q.addParamValue('minCellArea',250,@isnumeric);        %minimum cell area (smaller are erased)

%STEP C : prepare seeds for watershedding
q.addParamValue('GaussianFilter',5,@isnumeric);        %smoothing of the original image to find local minima within cells 
q.addParamValue('minDepth',5,@isnumeric);             %minimum accepted depth for a local minimum

%STEP E : treatment of long cells
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

%blubb start
% enforce conncetion
%propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');



%blubb end


propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');
[~,idx] = max([propsMaskImage.Area]);
%blubb A_cropMaskImage = propsMaskImage(idx).Image;                               %cropped mask
%bb = propsMaskImage(idx).BoundingBox;
%ROI_segmentation = [ceil(bb([2 1])) ceil(bb([2 1]))+bb([4 3])-[1 1]];   %ymin xmin ymax xmax
A_cropMaskImage=labelMaskImage;
%bb = [1 1 size(A_cropMaskImage,1) size(A_cropMaskImage,2)];
%increment a little bit, otherwise background calculation of fluor is going
%to fail
myIncr=15;
ROI_segmentation=[1+myIncr 1+myIncr size(A_cropMaskImage,1)-myIncr size(A_cropMaskImage,2)-myIncr];
A_cropMaskImage=A_cropMaskImage(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4));
%blubb

%
A_cropPhImage = O_PhImageFilt(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4)); %cropped ph image image

if q.Results.saveSteps
    savePNGofImage(A_maskImage,'A_maskImage',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP B : find edges 

B_negPh = imcomplement(A_cropPhImage);
se = strel('disk',1);
B_negPhErode = imerode(B_negPh, se);                                       %Morphological reconstruction ...
B_negPh = imreconstruct(B_negPhErode,B_negPh);                             %... in the mask of the original image
B_negPh = imdilate(B_negPh, se);                                           %Allows a better edge determination at contacts.
B_edgeImage1 = edge(B_negPh,'log',0,q.Results.LoG_Smoothing);                      %edge detection by smoothing (laplacian of gaussian ensures closed edges)

%suppress noisy surroundings
B_edgeImage2 = B_edgeImage1 & A_cropMaskImage;
B_fillEdgeImage2 = imfill(B_edgeImage2,'holes');
B_fillEdgeImage2 = bwareaopen(B_fillEdgeImage2,q.Results.minCellArea,4);           %suppress small stuff 
B_edgeImage2 = B_edgeImage1 & B_fillEdgeImage2;                            %keeps only edges that do not own to small stuff
B_edgeImage2 = bwareaopen(B_edgeImage2,30);                                %remove small loops related to intracell variations

if q.Results.saveSteps
savePNGofImage(B_edgeImage1,'B_edgeImage1',q.Results.saveDir);
savePNGofImage(B_edgeImage2,'B_edgeImage2',q.Results.saveDir);
savePNGofImage(B_fillEdgeImage2,'B_fillEdgeImage2',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP C : prepare seeds for watershedding

%makes a clean filled cells mask
C_smoothPh = imfilter(A_cropPhImage,fspecial('gaussian',q.Results.GaussianFilter,q.Results.GaussianFilter));
C_localMinPh = imextendedmin(C_smoothPh,q.Results.minDepth) & B_fillEdgeImage2;       %local minima in the mask      
C_cellMask = imfill(B_edgeImage2,find(C_localMinPh));                       
C_cellMask = bwmorph(C_cellMask,'open');                                    %clean cells mask

%shrinking steps to cut some cells
C_seeds1 = C_cellMask & ~B_edgeImage2;                                      %thins cells by removing the edge
C_seeds2 = bwmorph(C_seeds1,'open');                                        %already cuts some cells
C_seeds2 = bwmorph(C_seeds2,'thin',inf);                                    %skeletinization

if q.Results.saveSteps
savePNGofImage(C_cellMask & ~C_localMinPh,'C_Mask and minima',q.Results.saveDir);
savePNGofImage(C_seeds1,'C_seeds1',q.Results.saveDir);
savePNGofImage(C_seeds2,'C_seeds2',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP D : suppress branch points of the skeleton
brchpts = PN_FindBranchPoints(C_seeds2);
C_seeds2(logical(brchpts)) = 0;
%some cleaning
C_seeds2 = bwmorph(C_seeds2,'spur',3);
C_seeds2 = bwareaopen(C_seeds2,10,8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP E : cut long cells which neck is deeper than neckDepth
continueToCut = true;
icut = C_seeds2;
while continueToCut
    [cellsToRemove cutPoints] = PN_CutLongCells(icut,C_cellMask,q.Results.neckDepth);
    if max(max(cutPoints))==0
        continueToCut = false;
    else
        cutPoints = bwmorph(cutPoints,'dilate',2);
        C_seeds2(cutPoints) = false; %cuts the long cells on the seeds image
        icut(cutPoints) = false;
        icut= icut & ~cellsToRemove;
    end
end

if q.Results.saveSteps
savePNGofImage(C_seeds2,'C_seeds2',q.Results.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP Z : final segmentation by watershedding

%prepare mask for watershedding
Z_maskToFill = bwmorph(C_cellMask,'dilate');
Z_background = ~Z_maskToFill;                                               %background

%prepare final seeds
Z_seeds = bwmorph(C_seeds2,'spur',3);
Z_seeds = bwareaopen(Z_seeds,4,8);

%prepare seeded landscape for watershedding
Z_d1 = -bwdist(Z_background);                                               %distance transform  
Z_d1 = imimposemin(Z_d1, Z_background | Z_seeds);                           %impose background and seeds

%watershedding
Z_segmentedImage = watershed(Z_d1);                     
Z_segmentedImage(Z_background) = 0;
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
