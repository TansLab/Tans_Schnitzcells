function [A_cropPhImage Z_segmentedImage ROI_segmentation] = PN_segphase(imageToSegment,varargin)

% INPUTS:
%   imageToSegment:    phase contrast images (multiple Z slices)
%   varargin:      parameters for image treatment indicated below
%
% OUTPUTS:
%   A_cropPhImage:  phase contrast image cropped (default is median filtered input image)
%   Z_segmentedImage:    segmented image (smaller size)
%   ROI_segmentation:   corner coordinates of the crop region in the original image



%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%
%reads parameter values or initialize with default values not too bad for
%E. Coli colonies on acryl gel pad with 100X phase contrast images

q = inputParser;
q.addrequired = ('imageToSegment',@isnumeric)

%STEP A : finds a global mask and crop the image
q.addParamValue= ('rangeFiltSize',35,@isinteger);       %dectection of interesting parts of the image
q.addParamValue= ('maskMargin',20,@isinteger);          %additional margin of the mask : enlarge if cells missing on the edge

%STEP B : find edges
q.addParamValue= ('LoG_Smoothing',2,@isscalar);         %smoothing amplitude of the edge detection filter
q.addParamValue= ('minCellArea',250,@isinteger);        %minimum cell area (smaller are erased)

%STEP C : prepare seeds for watershedding
q.addParamValue= ('GaussianFilter',5,@isscalar);        %smoothing of the original image to find local minima within cells 
q.addParamValue= ('minDepth',5,@isinteger);             %minimum accepted depth for a local minimum

%STEP E: treatment of long cells
q.addParamValue= ('cutCellsWidth',3,@isscalar);         %minimum neck width to cut a too long cell

%saving images
q.addParamValue= ('saveDir',pwd,@isdir);                %directory where intermediate treatment images are saved
q.addParamValue= ('saveSteps',false,@islogical);            %indicate if you want to save intermediate images

q.parse(varargin{:});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP O (O the letter, not 0 the zero !): suppress shot noise
O_PhImageFilt = medfilt2(imageToSegment);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP A : find a global mask and crop the image

A_maskImage = rangefilt(O_PhImageFilt,true(rangeFiltSize));                 %detect zones of sufficient intensity variations
A_maskImage = im2bw(A_maskImage,graythresh(A_maskImage));                   %threshold to black and white
A_maskImage = imclose(A_maskImage,true(maskMargin));                        %enlarge mask         
labelMaskImage = bwlabel(A_maskImage);                                      %only keep the biggest connected part of the mask
propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');
[dumb idx] = max([propsMaskImage.Area]);
A_cropMaskImage = propsMaskImage(idx).Image;    %cropped mask
bb = propsMaskImage(idx).BoundingBox;
A_cropPhImage = imcrop(O_PhImageFilt, bb-[0 0 1 1]);   %cropped ph image image
ROI_segmentation = [floor(bb([1 2])) floor(bb([1 2]))+bb([3 4])-[1 1]];  

clear labelMaskImage propsMaskImage dumb idx
savePNGofImage(A_maskImage,'A_maskImage',saveDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP B : find edges 

B_negPh = imcomplement(A_cropPhImage);
se = strel('disk',1);
B_negPhErode = imerode(B_negPh, se);                                       %morphological reconstruction 
B_negPh = imreconstruct(B_negPhErode,B_negPh);                             %in the mask of the original image
B_negPh = imdilate(B_negPh, se);                                           %allows a better edge determination at contacts
B_edgeImage1 = edge(B_negPh,'log',0,LoG_Smoothing);                        %edge detection by smoothing (laplacian of gaussian ensures closed edges)

%suppress noisy surroundings
B_edgeImage2 = B_edgeImage1 & A_cropMaskImage;
B_fillEdgeImage2 = imfill(B_edgeImage2,'holes');
B_fillEdgeImage2 = bwareaopen(B_fillEdgeImage2,minCellArea,4);             %suppress small stuff 
B_edgeImage2 = B_edgeImage1 & B_fillEdgeImage2;                            %keeps only edges that do not own to small stuff
B_edgeImage2 = bwareaopen(B_edgeImage2,30);                                %remove small loops related to intracell variations

savePNGofImage(B_edgeImage1,'B_edgeImage1',saveDir);
savePNGofImage(B_edgeImage2,'B_edgeImage2',saveDir);
savePNGofImage(B_fillEdgeImage2,'B_fillEdgeImage2',saveDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP C : prepare seeds for watershedding

%makes a clean filled cells mask
h = fspecial('gaussian',GaussianFilter,GaussianFilter);
C_smoothPh = imfilter(A_cropPhImage,h);                                    %gaussian filtering of phimage
C_localMinPh = imextendedmin(C_smoothPh,minDepth) & B_fillEdgeImage2;      %local minima in the mask      
C_cellMask = imfill(B_edgeImage2,find(C_localMinPh));                      %clean cells mask

%shrinking steps to cut some cells
C_seeds1 = C_cellMask & ~B_edgeImage2;                                     %thins cells by removing the edge
C_seeds2 = bwmorph(C_seeds1,'open');                                       %already cuts some cells

savePNGofImage(C_cellMask & ~C_localMinPh,'Mask and minima',saveDir);
savePNGofImage(C_seeds1,'C_seeds1',saveDir);
savePNGofImage(C_seeds2,'C_seeds2',saveDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP D : treatment of kinky cells : to come.
%solidity is a first screen but more local criterion should be found


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP E: cut long cells with the threshold neck size cutCellsWidth

continue2cut = true; 
E_seedsCutShort = C_seeds2;
while continue2cut                                                          %cuts one connex domain per round, looped until no more cut found
    E_cutterImage = PN_CutLongCells(cutCellsWidth,E_seedsCutShort);
    E_cutterImage = bwmorph(E_cutterImage,'dilate',cutCellsWidth+2); 
    E_seedsCutShort(E_cutterImage) = false;
    continue2cut = ~isempty(find(E_cutterImage));
end

savePNGofImage(E_seedsCutShort,'E_seedsCutShort',saveDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP Z : final segmentation by watershedding

%prepare seeds for watershedding (it could be possible to begin from a new edge image) 
Z_seedsFinal = bwareaopen(E_seedsCutShort,minCellArea*0.8);                 %suppress small seeds which artificially cut within cell

%prepare mask for watershedding
Z_maskToFillClean = bwmorph(C_cellMask,'close');                            %prepare mask from which watershedding is performed
Z_maskToFillClean = bwareaopen(Z_maskToFillClean,minCellArea,4);
Z_maskToFill = bwmorph(Z_maskToFillClean,'dilate');                         %final mask 
Z_background = ~Z_maskToFill;                                               %background

%prepare landscape and seeds  for watershedding
Z_d1 = -bwdist(Z_background);                                               %distance transform  
Z_d1 = imimposemin(Z_d1, Z_background | Z_seedsFinal);                      %impose background and seeds

%watershedding
Z_w0 = watershed(Z_d1);                     
Z_w0(Z_background) = 0;
Z_w1 = bwareaopen(Z_w0,10);
Z_segmentedImage = bwlabel(Z_w1);                                           %segmentation is finished here

savePNGofImage(Z_maskToFillClean,'Z_maskToFillClean',saveDir);
savePNGofImage(Z_maskToFill & ~Z_seedsFinal,'Mask and seeds',saveDir);
savePNGofImage(Z_segmentedImage,'Z_segmentedImage',saveDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SEGMENTATION END   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
if max(max(Z_segmentedImage))==0                                           %let user know when no cells were found
  disp([' * !! WATCH OUT !! no cells found on this frame...']);
end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function savePNGofImage(image, name, saveDir)
if saveSteps
    filename = [saveDir,name,'.png'];
    if isempty(image)
        image = [ 0 ]; % will get error if image is empty
    end
    DJK_writeSegImage(image, filename);
    %disp(['          * written : ' name '.png']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save images of intermediate steps
function saveTIFFofPhaseImage(image, name, saveDir)
if saveSteps
    filename = [saveDir,name,'.tif'];
    if isempty(image)
        image = [ 0 ]; % will get error if image is empty
    end
    imwrite(uint16(image), [saveDir,name,'.tif'], 'TIFF', 'Compression', 'none');
    %disp(['          * written : ' name '.tif']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%