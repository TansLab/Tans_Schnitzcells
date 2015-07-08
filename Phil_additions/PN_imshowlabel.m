function outim = PN_imshowlabel(p,L,rect,Lp,rectp,varargin)
% function outim = PN_imshowlabel(p,L,rect,Lp,rectp,varargin)
%
% PN_imshowlabel is a modified version of DJK_imshowlabel which provides
% visual aid for segmentation correction ("suspicious cell detection")
% - too small cells are shown in white
% - cells with potential overlap with 2 cells in the former image are
% shown in white
%
% When rect, Lp and rectp are set to [], suspicious cell detection will not
% be activated. This should be done when running a "preliminary analysis".
% (Because this analysis relies on previous frames being present.)
% (In practice, susp. cell detection is activated if Lp is unequal to zero)
%
% OUTPUT
% 'outim'             color image
%
% REQUIRED ARGUMENTS:
% 'p'                   general movie info var
% 'L'                   seg image
%
% When the following vars are not set to [], the suspicious cell detection
% algorithm will be activated. 
% 'rect'                [MW: I'm suspecting this is the recteangle defining
%                       the location of the selected area from the ph image]
% 'Lp'                  [MW: I'm suspecting this is the L from previous
%                       frame]
% 'rectp'               [MW: I'm suspecting this is idem to rect, but from
%                       from previous frame.]
% %% MW 2014/08/29 - TODO: I'm not sure whether these descriptions are
% %%                 correct.
%
%
% OPTIONAL ARGUMENTS:
% 'phaseImage'        phase image of same size as seg, will be shown as
%                     background
% 'randomize' = 0     no randomizing of colormap (default:1)
% 'customColors'      provide your own colormap
% TODO: 'recalcCellNumbers' array for which region properties have to be
%                     recalculated because cell was updated. default: all
%                     cells. taking fewer cells will speed up process
% showNr              if set, prints cell idxs for this frame on cells.
%
%
%
% function outim = PN_imshowlabel(L,rect,Lp,rectp,varargin)

% elapsed time for frame 444 in 2012-05-08. 390 cells




%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------

% Settings
fractionbelowwhite = 0.2; % cells which are a fraction of <fractionbelowwhite>
                         % smaller than the median are marked white.
numRequiredArgs = 5; functionName = 'PN_imshowlabel'; p_internal = struct;

DISPLAY_REGION_ALPHA = 0.5; % alpha value for overlaying regions over phase img

% Input processing
if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)))
    errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
    error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
            error(errorMessage);
        end
        fieldName = DJK_schnitzfield(varargin{i});
        p_internal.(fieldName) = varargin{i+1};
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------

if existfield(p_internal,'phaseImage') & length(p_internal.phaseImage)>0
    addPhaseImage = true;
else
    addPhaseImage = false; % MW 
end

if ~existfield(p_internal,'randomize')
    p_internal.randomize = 1;
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

if islogical(L),
    L = double(L);
end;

%start time (calculate cost)
%tic

% If previous (preceeding) frame provided, perform suspicious cell. detection
if ~isempty(Lp)

    %Detect cells that have certainly been divided in the preceeding image
    [rw clm] = find(Lp);
    center_old = round([mean(clm) mean(rw)]);
    [rw clm] = find(L);
    center_new = round([mean(clm) mean(rw)]);
    pad_motion = center_new + [rect(2) rect(1)] - [rectp(2) rectp(1)] - center_old; 
    %Positions of the centroids of the former image in the new image
    propLp = regionprops(Lp,'Centroid');
    centroids = cat(1, propLp.Centroid);
    centroids = round(centroids + ones(size(propLp))*([rectp(2)-1 rectp(1)-1] + pad_motion));
    imcentroids = zeros(max(rect(3),rectp(3)),max(rect(4),rectp(4)));
    linearInd = sub2ind(size(imcentroids), centroids(:,2), centroids(:,1));
    imcentroids(linearInd) = 1;
    imcentroids = imcentroids(rect(1):rect(3),rect(2):rect(4));
    %find domains which have 2 or more former centroids
    Ltemp = L;
    Ltemp(logical(imcentroids)) = 0;
    propLtemp = regionprops(Ltemp,'EulerNumber');
    ideul = find([propLtemp.EulerNumber]<0); %label of cells which may have been re-merged
    %detect too small cells
    propL = regionprops(L,'Area');
    characSize = median([propL.Area]);
    idsmall = find([propL.Area] < characSize*fractionbelowwhite); %label of cells which may be too small 
                    % ^ fractionbelowwhite is the fraction of the avg under which the cell is colored white
    %create a logical of suspicious cells
    Lsuspicious = zeros(size(L));
    % If there are suspicous cells, mark them on the Lsuspicious image.
    allSuspiciousLabels = union(ideul,idsmall);
    if ~isempty(allSuspiciousLabels)
        for ii = allSuspiciousLabels
            Lsuspicious(L==ii)=1;
        end
    end

    % MW: make them stand out some more using checkerboard
    sqsize = 4; %square size
    mycheckerboard = (checkerboard(sqsize,ceil(size(Lsuspicious,1)/sqsize),ceil(size(Lsuspicious,2)/sqsize)) > 0.5);
    for i = 1:size(Lsuspicious,1)
        for j = 1:size(Lsuspicious,2)
            if mycheckerboard(i,j)==0
                Lsuspicious(i,j) = 0;
            end
        end
    end

    % MW (removed code, this statement can be removed)
    
end % note suspicious cell detection algorithm has a 2nd part below

% elapsed time: 0.33
%stop1=toc

%creation of a rgb colormap
% Use custom colormap if it is set
if isfield(p_internal,'customColors')
    mymap = p_internal.customColors;
    L2=uint8(L);
else
    % L2 has every non-background blob in range [2,256] and 
    % sets background to one, corresp. to first entry in mymap
    L2 = mod(L,255)+2;
    L2(L==0) = 1;
    % M is the maximum color table entry, at most 256 colors
    M = min(max2(L)+2,256);
    % create a color map
    mymap = DJK_hsv(M); % DJK 071207
    % explicitly set the colormap's first entry to black for background
    mymap(1,:)=[0 0 0];
    if p_internal.randomize
      % get sequence of random integers in range [1,maxcolors-1]
      [s,I] = sort(rand(M-1,1));  
      % randomly reorder mymap color entries [2,maxcolors]
      mymap(2:end,:) = mymap(I+1,:);
    end
    mymap = [mymap ; 1 1 1]; %add white
end

% 2nd part of suspicious cell detection
if ~isempty(Lp)
    % elapsed time: 0.34
    %stop2=toc
    L2(logical(Lsuspicious)) = M+1;
    L2(logical(imdilate(imcentroids,strel('square',2)))) = 0; % costs 0.017 sec    
end

Lrgb = ind2rgb(L2,mymap); % costs 0.04 sec

% elapsed time: 0.40
%stop3=toc
if existfield(p, 'showPerim') && p.showPerim % show cell outlines        
    % MW 2014/12        
    
    % Get phaseimg
    phaseImg = double(p_internal.phaseImage);
    phaseImg = DJK_scaleRange(phaseImg, [max(max(phaseImg)) min(min(phaseImg))], [0 1]);
    phaseImg = (1-phaseImg); % negative
    
    % Make 3d phaseimg
    outim=zeros([size(phaseImg),3]); % phase as base          
    outim(:,:,1) = phaseImg;
    outim(:,:,2) = phaseImg;
    outim(:,:,3) = phaseImg;  
    
    % Get perimeter image
    %perimImg = bwperim(L).*L; % make outline image; *.L colors the perims
    perimImg = bwperim(L);
    %perimImg = imdilate(perimImg,strel('disk',3,4)); % imdilate use TODO can be optimized
    perimImg = perimImg.*L;
    
    perimImg = ind2rgb(perimImg+1,mymap); % map indices to colours
          
    % Put perimeters in phase image
    nonZeroIdx = find(perimImg>0);
    outim(nonZeroIdx)=perimImg(nonZeroIdx);  

    % IF YOU WANT WHITE LINES:
    %{
    % Get phaseimg
    phaseImg = double(p_internal.phaseImage);
    phaseImg = DJK_scaleRange(phaseImg, [max(max(phaseImg)) min(min(phaseImg))], [0 1]);
    phaseImg = (1-phaseImg); % negative

    % Get perimeter
    perimImg = bwperim(L);%   
    
    % Create output image
    outim = phaseImg;    % original phase
    nonZeroIdx = find(perimImg>0); % get perimeter indexes
    outim(nonZeroIdx)=1; % mark them in phase img
    %}
    
    % Edit MW - adds green marker if frame is approved
    outim = MW_stampit(outim,p);
    
    if nargout == 0 && isempty(Lp)
        imshow(outim);
    end
    
elseif addPhaseImage % costs 0.045 sec
    %{
    rgb = 0.5 * Lrgb; % rgb = 0.5 * Lrgb; % MW here alpha set, TODO
    bwscreen = double(p_internal.phaseImage); % bwscreen = 0.5 * bwscreen / max(max(bwscreen));
    bwscreen = DJK_scaleRange(bwscreen, [max(max(bwscreen)) min(min(bwscreen))], [0 1]);
    bwscreen = DJK_scaleRange(bwscreen, [0.25 1], [0 0.5]);
    rgb(:,:,1) = rgb(:,:,1) + bwscreen;
    rgb(:,:,2) = rgb(:,:,2) + bwscreen;
    rgb(:,:,3) = rgb(:,:,3) + bwscreen;
    %}
    
    % (From http://en.wikipedia.org/wiki/Alpha_compositing)   
    outim = Lrgb.*DISPLAY_REGION_ALPHA; % region image as base    
    
    phaseImg = double(p_internal.phaseImage);
    phaseImg = DJK_scaleRange(phaseImg, [max(max(phaseImg)) min(min(phaseImg))], [0 1]);
    phaseImg = (1-phaseImg).*(1-DISPLAY_REGION_ALPHA); % (1-img) can make img negative
    
    outim(:,:,1) = outim(:,:,1) + phaseImg;
    outim(:,:,2) = outim(:,:,2) + phaseImg;
    outim(:,:,3) = outim(:,:,3) + phaseImg;    
    
    % elapsed time: 0.45
    %stop4a=toc
    
    % Edit MW - adds green marker if frame is approved
    outim = MW_stampit(outim,p); % TODO CHECK MW
    
    if nargout == 0 && isempty(Lp)
        imshow(outim);
    end
    
else
    
    outim = Lrgb;
    outim = MW_stampit(outim,p); % TODO CHECK MW

    if nargout == 0 && isempty(Lp)
        imshow(outim);
    end
    
    %stop4b=toc
end

if isfield(p,'showNr')
    FONTSIZE=8;
    halfFontSize=round(FONTSIZE/2);
    
    propL = regionprops(L,'Centroid');
    
    for i = 1:numel(propL)
        text(propL(i).Centroid(1)-halfFontSize,propL(i).Centroid(2)-halfFontSize,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[1,1,1],'FontWeight','bold')
    end
end

%elapsed time: 0.45
%stop5=toc
