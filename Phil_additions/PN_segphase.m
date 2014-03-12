function [A_cropPhImage, Z_segmentedImage, ROI_segmentation] = PN_segphase(imageToSegment,varargin)

% INPUTS:
%   imageToSegment:    greyscale image
%   varargin:      parameters for image treatment indicated below
%
% OUTPUTS:
%   A_cropPhImage:  phase contrast image cropped (default is median filtered input image)
%   Z_segmentedImage:    segmented image (smaller size)
%   ROI_segmentation:   corner coordinates of the crop region in the original image


%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%
%Reads parameter values or initializes with default values not too bad for
%E. Coli colonies on acryl gel pad with 100X phase contrast images.
%This is highly simplified by inputParser in more recent matlab releases.
q = struct; 

numRequiredArgs = 1;
if (nargin < 1) | (mod(nargin,2) == 0)
  errorMessage = sprintf ('%s\n%s\n%s\n','Error using PN_segphase:',...
      '    Invalid input arguments.','    Try "help PN_segphase".');
  error(errorMessage);
end


numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
          'Error using ==> PN_segphase_temp:',...
          '    Invalid property ', num2str(varargin{i}), ...
          ' is not (needs to be) a string.',...
          '    Try "help PN_segphase_temp".');
      error(errorMessage);
    end
    fieldName = varargin{i};
    q.(fieldName) = varargin{i+1};
  end
end

%%%%%%%%%%%default values

%STEP A : finds a global mask and crop the image
if ~existfield(q,'rangeFiltSize')                         %typical area for dectection of interesting features of the image
    q.rangeFiltSize = 35;       
end
if ~existfield(q,'maskMargin')                            %additional margin of the mask : enlarge if cells missing on the edge
    q.maskMargin = 20;       
end
%STEP B : find edges
if ~existfield(q,'LoG_Smoothing')                         %smoothing amplitude of the edge detection filter                        
    q.LoG_Smoothing = 2;       
end
if ~existfield(q,'minCellArea')                           %minimum cell area (objects smaller than that will be erased)                          
    q.minCellArea = 250;       
end
%STEP C : prepare seeds for watershedding
if ~existfield(q,'GaussianFilter')                        %smoothing of the original image to find local minima within cells                         
    q.GaussianFilter = 5;       
end
if ~existfield(q,'minDepth')                              %minimum accepted depth for a local minimum                         
    q.minDepth = 5;       
end
%STEP E: treatment of long cells
if ~existfield(q,'cutCellsWidth')                         %minimum neck width to cut a too long cell
    q.cutCellsWidth = 3.5;       
end
%saving images
if ~existfield(q,'saveSteps')                             %indicate if you want to save intermediate images
    q.saveSteps = false;    
end
if ~existfield(q,'saveDir') & q.saveSteps                 %directory where intermediate treatment images are saved
    error('Indicate a directory to save files')       
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP O (O the letter, not 0 the zero !): suppress shot noise
O_PhImageFilt = medfilt2(imageToSegment);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP A : find a global mask and crop the image

A_maskImage = rangefilt(O_PhImageFilt,true(q.rangeFiltSize));              %detect zones of sufficient intensity variations
A_maskImage = im2bw(A_maskImage,graythresh(A_maskImage));                  %threshold to black and white
A_maskImage = imclose(A_maskImage,true(q.maskMargin));                     %enlarge mask         
labelMaskImage = bwlabel(A_maskImage);                                     %only keep the biggest connected part of the mask
propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');
[dumb idx] = max([propsMaskImage.Area]);
A_cropMaskImage = propsMaskImage(idx).Image;                               %cropped mask
bb = propsMaskImage(idx).BoundingBox;
A_cropPhImage = imcrop(O_PhImageFilt, bb-[0 0 1 1]);                       %cropped ph image image
ROI_segmentation = [ floor(bb([2 1])) floor(bb([2 1]))+bb([4 3])-[1 1] ];   %ymin xmin ymax xmax

if q.saveSteps
    savePNGofImage(A_maskImage,'A_maskImage',q.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP B : find edges 

B_negPh = imcomplement(A_cropPhImage);
se = strel('disk',1);
B_negPhErode = imerode(B_negPh, se);                                       %Morphological reconstruction ...
B_negPh = imreconstruct(B_negPhErode,B_negPh);                             %... in the mask of the original image
B_negPh = imdilate(B_negPh, se);                                           %Allows a better edge determination at contacts.
B_edgeImage1 = edge(B_negPh,'log',0,q.LoG_Smoothing);                      %edge detection by smoothing (laplacian of gaussian ensures closed edges)

%suppress noisy surroundings
B_edgeImage2 = B_edgeImage1 & A_cropMaskImage;
B_fillEdgeImage2 = imfill(B_edgeImage2,'holes');
B_fillEdgeImage2 = bwareaopen(B_fillEdgeImage2,q.minCellArea,4);           %suppress small stuff 
B_edgeImage2 = B_edgeImage1 & B_fillEdgeImage2;                            %keeps only edges that do not own to small stuff
B_edgeImage2 = bwareaopen(B_edgeImage2,30);                                %remove small loops related to intracell variations

if q.saveSteps
savePNGofImage(B_edgeImage1,'B_edgeImage1',q.saveDir);
savePNGofImage(B_edgeImage2,'B_edgeImage2',q.saveDir);
savePNGofImage(B_fillEdgeImage2,'B_fillEdgeImage2',q.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP C : prepare seeds for watershedding

%makes a clean filled cells mask
h = fspecial('gaussian',q.GaussianFilter,q.GaussianFilter);
C_smoothPh = imfilter(A_cropPhImage,h);                                    %gaussian filtering of phimage
C_localMinPh = imextendedmin(C_smoothPh,q.minDepth) & B_fillEdgeImage2;    %local minima of greyscale within the mask      
C_cellMask = imfill(B_edgeImage2,find(C_localMinPh));                      %here we have a clean cells mask

%shrinking steps to cut some cells
C_seeds1 = C_cellMask & ~B_edgeImage2;                                     %to thin cells by removing the edge
C_seeds2 = bwmorph(C_seeds1,'open');                                       %already cuts some cells

if q.saveSteps
savePNGofImage(C_cellMask & ~C_localMinPh,'C_Mask and minima',q.saveDir);
savePNGofImage(C_seeds1,'C_seeds1',q.saveDir);
savePNGofImage(C_seeds2,'C_seeds2',q.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP D : treatment of kinky cells : to come.
%solidity is a first screen but more local criterion should be found


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP E: cut long cells with the threshold neck size cutCellsWidth

continue2cut = true; 
E_seedsCutShort = C_seeds2;
while continue2cut                                                         %cuts one connex domain per round, looped until no more cut found
    E_cutterImage = PN_CutLongCells(q.cutCellsWidth,E_seedsCutShort);
    E_cutterImage = bwmorph(E_cutterImage,'dilate',q.cutCellsWidth+2);     
    E_seedsCutShort(E_cutterImage) = false;                                %creates elements to cut the seeds
    continue2cut = ~isempty(find(E_cutterImage));
end

if q.saveSteps
savePNGofImage(E_seedsCutShort,'E_seedsCutShort',q.saveDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP Z : final segmentation by watershedding

%prepare seeds for watershedding (it could be possible to begin from a new edge image) 
Z_seedsFinal = bwareaopen(E_seedsCutShort,q.minCellArea*0.8);               %suppress small seeds which result in intracell cutting

%prepare mask for watershedding
Z_maskToFillClean = bwmorph(C_cellMask,'close');                            %prepare mask from which watershedding is performed
Z_maskToFillClean = bwareaopen(Z_maskToFillClean,q.minCellArea,4);
Z_maskToFill = bwmorph(Z_maskToFillClean,'dilate');                         %final mask 
Z_background = ~Z_maskToFill;                                               %background

%prepare landscape and seeds  for watershedding
Z_d1 = -bwdist(Z_background);                                               %distance transform of the final mask 
Z_d1 = imimposemin(Z_d1, Z_background | Z_seedsFinal);                      %impose background and seeds

%watershedding
Z_w0 = watershed(Z_d1);                     
Z_w0(Z_background) = 0;
Z_w1 = bwareaopen(Z_w0,40);
Z_segmentedImage = bwlabel(Z_w1);                                           %segmentation is finished here

if q.saveSteps
savePNGofImage(Z_maskToFillClean,'Z_maskToFillClean',q.saveDir);
savePNGofImage(Z_maskToFill & ~Z_seedsFinal,'Z_Mask and seeds',q.saveDir);
savePNGofImage(Z_segmentedImage,'Z_segmentedImage',q.saveDir);
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
    %disp(['          * written : ' name '.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
