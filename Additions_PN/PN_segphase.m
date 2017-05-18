% function [A_cropPhImage, Z_segmentedImage, ROI_segmentation] = PN_segphase(p, imageToSegment,varargin)
%
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
% 
% OPTIONAL INPUTS: (MW 2015/04)
%   useFullImage:   Default value 1; if set to 0 allows you to handle the
%                   image uncropped. (Mother function 
%                   PN_segmoviephase_3colors.m accepts p.useFullImage).
%   p.segmentationFigures 
%                   if this is a valid field figures with segmentation
%                   steps will be shown (independent from saving figures).
%
% TODO MW 2015/01
% Check whether hard-coded algorithm constants in this algorithm can be
% made into parameters with default values that can be adjusted.
%
% NOTE MW
% Also added "p" parameter as input; this is redundant with varargin/inputsOfSegmentation
% 
function [A_cropPhImage, Z_segmentedImage, ROI_segmentation] = PN_segphase(p,imageToSegment,varargin)


%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%
%%Reads parameter values or initializes with default values not too bad for
%E. Coli colonies on acryl gel pad with 100X phase contrast images.
%
% Set varargin={} to manually execute this part

q = inputParser;
q.addRequired('imageToSegment', @isnumeric);

%STEP A : finds a global mask and crop the image
q.addParamValue('rangeFiltSize',35,@isnumeric);       %dectection of interesting parts of the image
q.addParamValue('maskMargin',5,@isnumeric);          %additional margin of the mask : enlarge if cells missing on the edge
q.addParamValue('useFullImage',0,@isnumeric);         % use full image for segmentation of ROI mask, MW TODO @islogical?

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

%% Parse user values
q.parse(imageToSegment, varargin{:});

%%
%{
% Example settings
tempVarargin={'rangeFiltSize', 35}
q.parse(imageToSegment, tempVarargin{:});
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP O (O the letter, not 0 the zero !): suppress shot noise
O_PhImageFilt = medfilt2(imageToSegment);

if isfield(p,'segmentationFigures')
    figure(1); imshow(O_PhImageFilt,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP A : find a global mask and crop the image


% Find a global mask if desired (see useFullImage param) and crop the image
A_maskImage = rangefilt(O_PhImageFilt,true(q.Results.rangeFiltSize));              %detect zones of sufficient intensity variations
A_maskImage = im2bw(A_maskImage,graythresh(A_maskImage));                  %threshold to black and white
A_maskImage=imfill(A_maskImage,'holes');                            % fill interior holes
A_maskImage = imdilate(A_maskImage,strel('disk',q.Results.maskMargin));                     %enlarge mask         
% ** NW2014-07: used to be 'imclose' instead of 'imdilate' but this didn't really enlarge the mask, only close it

labelMaskImage = logical(A_maskImage);                                     %only keep the biggest connected part of the mask



% determine region of interest
if isfield(p,'customColonyCenter') % MW addition
    
    propsMaskImage = regionprops(labelMaskImage,'BoundingBox','Image','Centroid');
    
    mycenter = p.customColonyCenter;
    
    distancestocenter = cellfun(@(x) sqrt((x(1)-mycenter(1)).^2 + (x(2)-mycenter(2)).^2), ...
        {propsMaskImage.Centroid});
   
    [~,idx] = min(distancestocenter);
    
else
    propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image','Centroid');
    [~,idx] = max([propsMaskImage.Area]);
end

if q.Results.useFullImage~=1  % default: crop image to ROI
    A_cropMaskImage = propsMaskImage(idx).Image;                               %cropped mask
    bb = propsMaskImage(idx).BoundingBox;
    ROI_segmentation = [ceil(bb([2 1])) ceil(bb([2 1]))+bb([4 3])-[1 1]];   %ymin xmin ymax xmax
else  % exception: skip cropping
    
    A_cropMaskImage=ones(size(imageToSegment));
    ROI_segmentation = [1 1 size(imageToSegment)];
    
    %crop a little bit, otherwise background calculation of fluor is going
    %to fail
    %{
    myIncr=15;
    ROI_segmentation=[1+myIncr 1+myIncr size(A_cropMaskImage,1)-myIncr size(A_cropMaskImage,2)-myIncr];
    A_cropMaskImage=A_cropMaskImage(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4));
    %}
end

A_cropPhImage = O_PhImageFilt(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4)); %cropped ph image image

if q.Results.saveSteps
    savePNGofImage(A_maskImage,'A_maskImage',q.Results.saveDir);
    %saveas(A_maskImage,'A_maskImage','png')
end

if isfield(p,'segmentationFigures')
    figure(1); imshow(A_cropPhImage,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP B : find edges 

B_negPh = imcomplement(A_cropPhImage);
se = strel('disk',1);
B_negPhErode = imerode(B_negPh, se);                                       %Morphological reconstruction ...
B_negPh = imreconstruct(B_negPhErode,B_negPh);                             %... in the mask of the original image
B_negPh = imdilate(B_negPh, se);                                           %Allows a better edge determination at contacts.
B_edgeImage1 = edge(B_negPh,'log',0,q.Results.LoG_Smoothing);              %edge detection by smoothing (laplacian of gaussian ensures closed edges)

%suppress noisy surroundings
B_edgeImage2 = B_edgeImage1 & A_cropMaskImage;
B_fillEdgeImage2 = imfill(B_edgeImage2,'holes');
B_fillEdgeImage2 = bwareaopen(B_fillEdgeImage2,q.Results.minCellArea,4);   %suppress small stuff 
B_edgeImage2 = B_edgeImage1 & B_fillEdgeImage2;                            %keeps only edges that do not own to small stuff
B_edgeImage2 = bwareaopen(B_edgeImage2,30);                                %remove small loops related to intracell variations


DE_boolean=0; % TODO make this cleaner
if DE_boolean
    %DE 2013-07-09. Issue with bright empty areas between the cells detected as
    % jagged regions. Annoying while manual segmentaion.
    % Make an additional mask and use it later to remove unwanted jagged regions.
    
    % threshold: if int<mean(int): int=0; this eliminates the most of the unwanted jagged regions
    B_negPh(B_negPh<mean(B_negPh(:)))=0;
    %make a binary mask from it:
    mask_fromintensity=(logical(B_negPh));
    %erode the mask, otherwise it's too tight and can cut the edges of the cells drastically (lysis :))
    se5 = strel('disk',10);
    mask_fromintensity2=imcomplement(imerode(imcomplement(mask_fromintensity),se5));
    
    % apply the mask obtained above to cut away the unwanted
    % regions from the old (B_edgeImage2) edge matrix, containing the edges of 
    % unwanted regions too:
    B_edgeImage3 = B_edgeImage2 .* mask_fromintensity2;
    
    % one can actually check the difference:
    % edge_difference=B_edgeImage2-B_edgeImage3; imagesc(edge_difference);
    % the rest goes unchanged:
    B_fillEdgeImage3 = imfill(B_edgeImage3,'holes');
    B_fillEdgeImage3 = bwareaopen(B_fillEdgeImage3,q.Results.minCellArea,4);  
    B_edgeImage3 = B_edgeImage3 & B_fillEdgeImage3;                            %keeps only edges that do not own to small stuff
    B_edgeImage3 = bwareaopen(B_edgeImage3,30);
    
    % IMPORTANT!!!! HERE I REPLACE B_edgeImage2 with B_edgeImage2 !!!!!
    B_edgeImage2 = B_edgeImage3;
end


if q.Results.saveSteps
savePNGofImage(B_edgeImage1,'B_edgeImage1',q.Results.saveDir);
savePNGofImage(B_edgeImage2,'B_edgeImage2',q.Results.saveDir);
savePNGofImage(B_fillEdgeImage2,'B_fillEdgeImage2',q.Results.saveDir);
end

if isfield(p,'segmentationFigures')
    figure(1); imshow(B_edgeImage1,[]);
    figure(2); imshow(B_edgeImage2,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP C : prepare seeds for watershedding

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

if isfield(p,'segmentationFigures')
    figure(1); imshow(C_cellMask,[]);
    figure(2); imshow(C_seeds2,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP D : suppress branch points of the skeleton

brchpts = PN_FindBranchPoints(C_seeds2);
C_seeds2(logical(brchpts)) = 0;
%some cleaning
C_seeds2 = bwmorph(C_seeds2,'spur',3);
C_seeds2 = bwareaopen(C_seeds2,10,8);

if isfield(p,'segmentationFigures')
    figure(1); imshow(C_seeds2,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP E : cut long cells which neck is deeper than neckDepth

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

if isfield(p,'segmentationFigures')
    figure(1); imshow(C_seeds2,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP Z : final segmentation by watershedding

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

% addition MW 2015/04, re-surpressing small cells (since after all
% operations, new artifacts might have been introduced).
% So note: this is the second time this operation is performed
Z_segmentedImage = MW_removesmallsegmented(Z_segmentedImage,q.Results.minCellArea); 

if q.Results.saveSteps
savePNGofImage(Z_maskToFill & ~Z_seeds,'Z_Mask and seeds',q.Results.saveDir);
savePNGofImage(Z_segmentedImage,'Z_segmentedImage',q.Results.saveDir);
end

if isfield(p,'segmentationFigures')
    figure(1); imshow(Z_segmentedImage,[]);
    p.showPhaseImage=0;
    figure(2); PN_imshowlabel(p,Z_segmentedImage,[],[],[]);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION END   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
if max(max(Z_segmentedImage))==0                                           %let user know when no cells were found
  disp([' * !! WATCH OUT !! no cells found on this frame...']);
end;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function savePNGofImage(image, name, saveDirectory)
    filename = [saveDirectory,name,'.png'];
    if isempty(image)
        image = [ 0 ]; % will get error if image is empty
    end
    DJK_writeSegImage(image, filename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
