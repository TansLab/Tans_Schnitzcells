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
% p.showNr            if set, prints cell idxs for this frame on cells.
% problemCells        (requires p.currentFrame to be set to current frame
%                     index) If p.problemCells is given, cells from that
%                     array will be highlighted (too). 
%                     problemcells = [schnitznr, framenr, labelnr; ..]
% p.showPerim         if this is a valid field, sohws outlines of cells
% p.showPhaseImage    if this field exists and is false, the phase image is
%                     hidden (per default phase image is shown and 
%                     p.showPhaseImage is set to 1).q
% p.slookup           if this field is set to a frame label to schnitz nr
%                     lookup table, schnitz nrs will be displayed.
% p.showLargeCentroids
%                     if this is a valid field the centroids from cells in
%                     the previous frame will be shown as large circles
%
%
%
% function outim = PN_imshowlabel(L,rect,Lp,rectp,varargin)

% elapsed time for frame 444 in 2012-05-08. 390 cells




%%-------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------

% Settings
fractionbelowwhite = 0.2; % cells which are a fraction of <fractionbelowwhite>
                         % smaller than the median are marked white.
numRequiredArgs = 5; functionName = 'PN_imshowlabel'; p_internal = struct;

DISPLAY_REGION_ALPHA = 0.5; % alpha value for overlaying regions over phase img

if ~isempty(Lp)
    assistedCorrection = 1;
else
    assistedCorrection = 0;
end

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

% If schnitznrs known, might as well map framelbls to schnitznrs
useSchnitzColors=0;
if isfield(p,'showNr') && isfield(p,'slookup')
    if (p.showNr==2)
        useSchnitzColors=1;
    end
end
%--------------------------------------------------------------------------


%%-------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------

if ~existfield(p,'showPhaseImage') 
    p.showPhaseImage = true;
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

sizeL = size(L);

%start time (calculate cost)
%tic

%% If previous (preceeding) frame provided, perform suspicious cell. detection
if assistedCorrection

    %Detect cells that have certainly been divided in the preceeding image
    [rw clm] = find(Lp);
    center_old = round([mean(clm) mean(rw)]);
    [rw clm] = find(L);
    center_new = round([mean(clm) mean(rw)]);
    pad_motion = center_new + [rect(2) rect(1)] - [rectp(2) rectp(1)] - center_old; 
    
    %Positions of the centroids of the former image in the new image
    propLp = regionprops(Lp,'Centroid');
    centroidsraw = cat(1, propLp.Centroid);
    centroids = round(centroidsraw + ones(size(propLp))*([rectp(2)-1 rectp(1)-1] + pad_motion));
    imcentroids = zeros(max(rect(3),rectp(3)),max(rect(4),rectp(4)));
    
    try
        linearInd = sub2ind(size(imcentroids), centroids(:,2), centroids(:,1));
    catch
        warning('Converting centroids to linear indices failed.');
        linearInd = NaN;
    end
    
    if ~isnan(linearInd)
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
    
        % Labels of suspicious cells
        allSuspiciousLabels = union(ideul,idsmall);
                
    else
        allSuspiciousLabels = [];
    end
end % more "if assistedCorrection" statements below

% also highlight problemcells, if desired.
if isfield(p,'problemCells');
    problemCellHighlighting = 1;
    
    if ~exist('allSuspiciousLabels','var')
        allSuspiciousLabels = [];
    end
    
    % If problemcells are known, we can label them in this image
    if ~isempty(p.problemCells) & isfield(p,'currentFrame')
        schnitzNrsProblemCells = find(p.problemCells(:,2)==p.currentFrame);
        allSuspiciousLabels = [ allSuspiciousLabels ...
                    p.problemCells(schnitzNrsProblemCells,3)'];    
    end
else 
    problemCellHighlighting = 0;
end

%% if either suspicious cell detection || problemcell highlighting,
% start highlighting cells.
if assistedCorrection || (problemCellHighlighting && (~isempty(allSuspiciousLabels)))        
    %create a logical of suspicious cells
    Lsuspicious = zeros(size(L));
        
    % If there are suspicous cells, mark them on the Lsuspicious image.
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

%% creation of a rgb colormap
% Use custom colormap if it is set
if useSchnitzColors
    mymap = p.customColors;
    L2=L;
    highestSchnitzIndx = max(p.slookup(size(p.slookup,1),:));
    highestCellFrameIndx = size(p.slookup,2); % = highest index in L
        % highestCellFrameIndx == size(p.customColors,1)-2
else
    % L2 has every non-background blob in range [2,256] and 
    % sets background to one, corresp. to first entry in mymap
    L2 = mod(L,255)+2;
    L2(L==0) = 1;
    % M is the maximum color table entry, at most 256 colors
    %highestCellFrameIndx = min(max2(L)+2,256); % white-color bug
    highestCellFrameIndx = max2(L)+2;
    % create a color map
    mymap = DJK_hsv(highestCellFrameIndx); % DJK 071207
    % explicitly set the colormap's first entry to black for background
    mymap(1,:)=[0 0 0];
    if p_internal.randomize
      % get sequence of random integers in range [1,maxcolors-1]
      [s,I] = sort(rand(highestCellFrameIndx-1,1));  
      % randomly reorder mymap color entries [2,maxcolors]
      mymap(2:end,:) = mymap(I+1,:);
    end
    mymap = [mymap; 1 1 1]; %duplicate map and add white
end

if useSchnitzColors
    % Create conversiontable conversionTable(framecellnr+1)=schnitznr
    lookupForThisFrame = p.slookup(p.currentFrame, :);
    
    % cells that are not mapped should become white
    lookupForThisFrame(lookupForThisFrame==0) = highestCellFrameIndx+1; 
    
    conversionTable = [0 lookupForThisFrame highestSchnitzIndx+1]; % add 0 for entry 0
        % adding extra index at end, which will correspond to color white later
    
    % Conversion table is applied later
    
    % Alternative way
    % possible but takes to long because contains loop
    % L = changem(L, conversionTable, 1:numel(conversionTable));    

end

%% 2nd part of suspicious cell detection || problemcellhighlighting
if assistedCorrection || (problemCellHighlighting && ~isempty(allSuspiciousLabels))
    % elapsed time: 0.34
    %stop2=toc
    L2(logical(Lsuspicious)) = highestCellFrameIndx+1;
    if assistedCorrection % only for susp. detection.
        L2(logical(imdilate(imcentroids,strel('square',2)))) = 0; % costs 0.017 sec    
    end
end


if useSchnitzColors
    
    % Debug code
    %{
    warning('Debug code activated.');
    sizeTable = size(conversionTable)
    maxL2 = max(L2(:)+1)
    Lindices = unique(L2(:))
        
    disp('Colors I''ll use:');
    colorIdx = conversionTable(Lindices+1)
    colors = mymap(colorIdx+1,:)+1
    
    % End debug code
    %}
    
    schnitzSegmentationMatrix = conversionTable(L2+1)+1;
    
    Lrgb = ind2rgb(schnitzSegmentationMatrix,mymap); 
        % 2x +1 since L contains 0 for blackground, and accordingly all
        % both conversionTable and colormap are shifted by 1.
        % Note that there's also the color white reserved at the end of the
        % colormap.
    
else 
    Lrgb = ind2rgb(L2,mymap); % costs 0.04 sec
end
    

% elapsed time: 0.40
%stop3=toc
if existfield(p, 'showPerim') && p.showPerim % show cell outlines        
    %% MW 2014/12                
    
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
    %perimImg = bwperim(L);
        % note that by definition touching stuff isn't recognized as
        % boundary
    % make lines thicker
    %perimImg = imdilate(perimImg,strel('disk',1,4)); % imdilate use TODO can be optimized            
        
    perimImg = L;
    
    cellBodies=zeros(size(L));
    
    for cellIdx = 1:max(max(L))
        
        workImg=zeros(size(L));
        workImg(L==cellIdx)=1;
        
        cellBodies = cellBodies+imerode(workImg,strel('disk',1,4));
    end
    
    % actually put in the lines by removing from original color-img areas
    % that are no line
    %perimImg = perimImg.*L;    
    perimImg(cellBodies>0) = 0;
    
    if useSchnitzColors
        perimImg = conversionTable(perimImg+1)+1;
        
        %perimImg = conversionTable(perimImg+1)+1;    
        perimImg = ind2rgb(perimImg,mymap); 
    else
        perimImg = ind2rgb(perimImg+1,mymap); % map indices to colours
    end
    
          
    % Put perimeters in phase image
    nonZeroIdx = perimImg(:,:,1)>0 | perimImg(:,:,2)>0 | perimImg(:,:,3)>0; % some RGB values have a 0 in them, e.g. red = [1,0,0]
    nonZeroIdx = [nonZeroIdx,nonZeroIdx,nonZeroIdx];
    
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
        %clf; % to prevent memory filling up
        imshow(outim);
        % hard print frame nr
        if isfield(p,'currentFrame'), text(17, 17, ['fr #' num2str(p.currentFrame)], 'FontSize', 12, 'Color', 'g'); end
    end
    
elseif p.showPhaseImage % costs 0.045 sec
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
        clf(gca); % to prevent memory filling up
        imshow(outim);
        % hard print frame nr
        if isfield(p,'currentFrame'), text(17, 17, ['fr #' num2str(p.currentFrame)], 'FontSize', 12, 'Color', 'g'); end
    end
    
else
    
    outim = Lrgb;
    outim = MW_stampit(outim,p); % TODO CHECK MW

    if nargout == 0 && isempty(Lp)
        clf(gca); % to prevent memory filling up
        imshow(outim);
        % hard print frame nr
        if isfield(p,'currentFrame'), text(17, 17, ['fr #' num2str(p.currentFrame)], 'FontSize', 12, 'Color', 'g'); end
    end
    
    %stop4b=toc
end

% For assisted correction, outim is resized
if assistedCorrection
    % outim = imresize_old(outim,p.res); % MW REMOVED LINE REMOVE THIS
    
    %imshow(outim,'InitialMagnification','fit');    % MW REMOVED LINE REMOVE THIS
    clf(gca); % to prevent memory filling up
    imshow(outim);

    % hard print frame nr
    %if isfield(p,'currentFrame'), text(40*p.res, 17, ['fr #' num2str(p.currentFrame)], 'FontSize', 12, 'Color', 'g'); end
    if isfield(p,'currentFrame'), text(17, 17, ['fr #' num2str(p.currentFrame)], 'FontSize', 12, 'Color', 'g'); end
end

if isfield(p,'showNr')
if p.showNr~=0
    FONTSIZE=8;
    halfFontSize=round(FONTSIZE/2);
    
    propL = regionprops(L,'Centroid');
    
    % print schnitznrs if lookup table available and option chosen
    if isfield(p,'slookup') && (p.showNr==2)
                
        for i = 1:numel(propL)
            % if lookup is available (might not be the case if manual 
            % correction to seg made)
            if ((size(p.slookup,2))>=i) % MW TODO
                % lookup nrs and print 
                schnitzNr = p.slookup(p.currentFrame,i);
                textx=propL(i).Centroid(1)-halfFontSize;
                texty=propL(i).Centroid(2)-halfFontSize;
                text(textx,texty,sprintf('%03d', schnitzNr),'FontSize',FONTSIZE,'Color',[1,1,1],'FontWeight','bold')                
            else
                % Tell users that other numbers might be incorrect too..
                warning(['Lookup of schnitznrs might be wrong, probably bc you corrected seg (attempted f= ' num2str(p.currentFrame) ' i=' num2str(i) ')..']);
            end
        end
    else
        % otherwise just print label
        for i = 1:numel(propL)
            textx=propL(i).Centroid(1)-halfFontSize;
            texty=propL(i).Centroid(2)-halfFontSize;
            text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[1,1,1],'FontWeight','bold')            
        end
    end
end

end

% Option below shows centroids from previous frame. The colonies tend to
% shift a bit because the field of view tends to shift a bit. The white
% circle is corrected for this, the grey one is not.
% This option is an addition to the black dots already shown, which are
% equivalent to the white circles.
if isfield(p,'showLargeCentroids')
    if p.showLargeCentroids & exist('centroidsraw','var')
            hold on; plot(centroidsraw(:,1),centroidsraw(:,2),'o','Color',[.7 .7 .7]);
            % pad_motion
            hold on; plot(centroidsraw(:,1) + pad_motion(1), centroidsraw(:,2) + pad_motion(2),'o','Color',[1 1 1]);
    end
end

if isfield(p,'showmubar')    
            
            imSize=size(phaseImg);
                        
            lengthBar = 1./p.micronsPerPixel;
            
            posx=[10,10+lengthBar];
            posy=[imSize(1)-10,imSize(1)-10];
            posytext=imSize(1)-27;
            
            %plot([10,10+lengthBar],[imSize(1)-10,imSize(1)-10],'LineWidth',5,'Color',[1 1 1]);
            hold on; plot(posx,posy,'LineWidth',5,'Color',[1 1 1]);
            text(posx(1),posytext,'1 µm','FontSize', 12,'Color',[1,1,1])
    
end

%elapsed time: 0.45
%stop5=toc
