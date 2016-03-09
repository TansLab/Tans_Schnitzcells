
% TODO MW 2015/01
% Investigate whether this function can be merged into PN_segphase as a
% functionality. Idem for NW_segphase_diffStrains.

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
q.addParamValue('maskMargin',5,@isnumeric);          %additional margin of the mask : enlarge if cells missing on the edge
q.addParamValue('useFullImage',0,@isnumeric);         % use full image for segmentation of ROI mask

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
A_maskImage=imfill(A_maskImage,'holes');                            % fill interior holes
A_maskImage = imdilate(A_maskImage,strel('disk',q.Results.maskMargin));                     %enlarge mask         
% ** NW2014-07: used to be 'imclose' instead of 'imdilate' but this didn't really enlarge the mask, only close it

labelMaskImage = logical(A_maskImage);                                     %only keep the biggest connected part of the mask
propsMaskImage = regionprops(labelMaskImage,'Area','BoundingBox','Image');
[~,idx] = max([propsMaskImage.Area]);
if q.Results.useFullImage~=1  % default: crop image to ROI
    A_cropMaskImage = propsMaskImage(idx).Image;                               %cropped mask
    bb = propsMaskImage(idx).BoundingBox;
    ROI_segmentation = [ceil(bb([2 1])) ceil(bb([2 1]))+bb([4 3])-[1 1]];   %ymin xmin ymax xmax
else  % exception: skip cropping
   % A_cropMaskImage=labelMaskImage;
   % %crop a little bit, otherwise background calculation of fluor is going
   % %to fail
   % myIncr=15;
   % ROI_segmentation=[1+myIncr 1+myIncr size(A_cropMaskImage,1)-myIncr size(A_cropMaskImage,2)-myIncr];
   % A_cropMaskImage=A_cropMaskImage(ROI_segmentation(1):ROI_segmentation(3),ROI_segmentation(2):ROI_segmentation(4));
   error('Rich medium & useFullImage will probably fail (extensive calculation time). Use default medium.')
end

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

%icut(logical(brchptsicut)) = 0; % removed NW2014-03 rich medium

%some cleaning
icut = bwmorph(icut,'spur',3);
icut = bwareaopen(icut,10,8);

while continueToCut
    [cellsToRemove cutPoints] = NW_CutLongCells_richMed(icut,C_cellMaskSmall,q.Results.neckDepth);

    % for debugging:
    %figure
    %imagesc(C_cellMaskSmall+2*imdilate(cutPoints,strel('disk',5)))
    
    if max(max(cutPoints))==0
        continueToCut = false;
    else
        % ********************************
        %create a divsion line -> similar to seg correction (updated NW
        %2014-03)
        
        %get coordinates of cutPoints (so far stored as points in a matrix)
        [idxrow, idxcol]=find(cutPoints>0);
        
        %label areas preliminarily to find which cell is to be cut
        C_prelimLabel=bwlabel(C_cellMaskSmall);
        
        for run=1:length(idxrow)          

            % for simplicity: convert notation to same as in PN_manual_kant
            cutx = idxrow(run); % CAREFUL! x stand here for rows!!!!
            cuty = idxcol(run);

            % extract cell that will be cut
            chosencolor = C_prelimLabel(cutx,cuty);
            cell = zeros(size(C_prelimLabel));
            cell(C_prelimLabel==C_prelimLabel(cutx,cuty)) = 1;%Lout(Lout==chosencolor);
            [fx,fy] = find(cell);
            xmin = max(min(fx)-5,1);
            xmax = min(max(fx)+5,size(cell,1));
            ymin = max(min(fy)-5,1);
            ymax = min(max(fy)+5,size(cell,2));
            subcell = cell(xmin:xmax, ymin:ymax); % subcell is only cell that will be cut

            % perim is perimeter of dilated cell
            perim = bwperim(imdilate(subcell,strel('disk',1)));
            % starting from cutPoitns-point, will increase a box until 2 sides are
            % found: this will be perims
            perims = zeros(size(perim));
            radp = 1;
            while max2(perims)<2 & radp<41
                 pxmin = max(cutx-xmin+1-radp,1);
                 pxmax = min(cutx-xmin+1+radp,size(perims,1));
                 pymin = max(cuty-ymin+1-radp,1);
                 pymax = min(cuty-ymin+1+radp,size(perims,2));
                 perims(pxmin:pxmax,pymin:pymax) = bwlabel(perim(pxmin:pxmax,pymin:pymax));
                 radp = radp+1;
            end
            % if indeed 2 sides are found, will cut
            if max2(perims)>1
                % kim is image with only clicked point drawn
                kim=zeros(size(subcell));
                kim(cutx-xmin+1,cuty-ymin+1)=1;

                % look for start of drawline
                kim1=kim;
                % increase size of kim untill it hits perims
                while ~any(any(kim1 & perims))
                    kim1=imdilate(kim1,strel('disk',1));
                end
                % randomly select first point as start of drawline
                [cut1x,cut1y]=find(kim1 & perims);

                % now go for end of drawline, first remove points of side of start from perims
                 color1=perims(cut1x(1),cut1y(1));
                 perims(perims==color1)=0;
                 kim2=kim;
                 while ~any(any(kim2 & perims))
                    kim2=imdilate(kim2,strel('disk',1));
                 end
                 % randomly select first point as end of drawline
                 [cut2x,cut2y]=find(kim2 & perims);
                 color2=perims(cut2x(1),cut2y(1));

                 % cut cell by drawing a thick(!) line (different to seg
                 % correction)
                 linethick=zeros(size(subcell));
                 linethick=drawline(linethick,[cut1x(1) cut1y(1)],[cut2x(1) cut2y(1)],1);
                 % dilate (thicken) line vertical to its major axis (the
                 % intuitive "thickening")
                 lineangle=regionprops(linethick,'Orientation'); % major axis orientation
                 lineangle=lineangle.Orientation;
                 SE=strel('line',4,lineangle+90); %structuring element is a line vertical to the division line
                 % a bit of an awkward way to thicken the strel line by 1
                 % (otherwise artefacts may appear since the thickened line
                 % can have holes)
                 mystrel=getnhood(SE);
                 mystrel=imdilate(mystrel,strel('disk',1));
                 SE=strel('arbitrary',mystrel);
                 % thicken division line
                 linethick=imdilate(linethick,SE);
                 subcell(linethick==1)=0;
                 % cutcell is original cell, but now with seperate cells different colors
                 cutcell = cell;
                 cutcell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4); % 4-connected

                 % remove original cell in Lout
                 C_prelimLabel(C_prelimLabel==C_prelimLabel(cutx,cuty))=0;

                 % first cell gets original color
                 C_prelimLabel(cutcell==1) = chosencolor;

                 % new cells get new color
                 for k = 2:max2(cutcell),
                     C_prelimLabel(cutcell==k) = max2(C_prelimLabel)+k-1;
                 end
            end

        % **********************************************
        
        end
        C_cellMaskSmall=C_prelimLabel>0;
        cutPointsDilated = bwmorph(cutPoints,'dilate',2);   %  is also original way how to cut cells
        icut(cutPointsDilated) = false;

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
