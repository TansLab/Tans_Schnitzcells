function [p, allLengthsOfBacteriaInPixels, allLengthsOfBacteriaInMicrons,allLengthsOfBacteriaInPixelsMW, allLengthsOfBacteriaInMicronsMW] = NDL_lengthforfillamentedcellsMW(p, frameRange) 
%%function [p,schnitzcells] = NDL_lengthforfillamentedcells(p) 
% 
% Function written by Nick de Lange and Martijn Wehrens.
% 2016/04
% 
% This function calculates length for bacteria using skeletons. This is
% useful when bacteria are very "curly" and fitting a curve to their shape
% does not yield succesfull results.
%
% Additionally, if p.generateStraigthenedBacteria is set, also images of
% straigthened bacteria will be generated using the function
% MW_straightenbacteria. (Note that this function outputs two parameters, the
% straigthened bacteria, and the length progression along the skeleton,
% since the transformation does not preserve equidistance of pixels 
% within the bacteria.)
%
% Input arguments:
% p.extraOutput         if this parameter is set, extra output will be 
%                       shown in plots and in command window.
% p.saveSummaryPlots    optional parameter that lets you decide to generate
%                       output on how skeleton looks to plots
%
% Output arguments:
%     % lengths
%     allLengthsOfBacteriaInPixels{framenr}(cellnum)
%     allLengthsOfBacteriaInMicrons{framenr}(cellnum)
%     
%     % Skeleton, area, additional data
%     allPixelAreaOfBacterium{framenr}(cellnum) 
%     allSkeletonXYpoleToPole{framenr}(cellnum) 
%     allMinX{framenr}(cellnum)
%     allMinY{framenr}(cellnum) 
%     allarray2{framenr}{cellnum}
%     alldistanceAlongSkeleton{framenr}{cellnum}
%     allextrapolatedDistanceEnds{framenr}(:,cellnum)
%
% Test easily with following command:
% >> NDL_lengthforfillamentedcells(p, ourSettings.frameRangeFull) 
%
%

% function parameters set by user:
AVERAGEBACTERIAWIDTH = 0.5; % In micron
DIRECTIONVALUECORRECT = 1.25; % Amount how much gets extrapolated in the 'correct' direction (times EXTRAPOLATIONLENGTH)
DIRECTIONVALUEINCORRECT = 0.4; % Amount how much gets extrapolated in the 'incorrect' direction (times EXTRAPOLATIONLENGTH)
% EXTRAPOLATIONLENGTH = 30; % In pixels - Dependant of pixel size --> CHANGE WHEN PIXEL SIZE IS DIFFERENT

% Error threshold values set by user:
ERRORVALUESKELETONLENGTH = 1; % In pixels - Gives error when extrapolation went non-optimal for skeletons larger than this number
ERRORFACTOREXTRAPOLATIONLENGTH = 1.2; % Gives error when extrapolation is at least this amount of times larger than EXTRAPOLATIONLENGTH
ERRORINTERSECT = 0.2; % Gives an error when the closest intersection between extrapolation and edge is still apart this number times EXTRAPOLATIONLENGTH
ERRORINTERSECTWARNING = 0.1; % Same as ERRORINTERSECT, but with lower threshold which can mean extrapolation is not optimal

% For smoothing of skeleton
SMOOTHELEMENTS = 8;

%frameRange = 241%unique([schnitzcells(:).frame_nrs]);
if exist('ourSettings', 'var') == 1
    frameRange = ourSettings.frameRangeFull; % Sets framerange to the full framerange provided in the Excel file
    warning('Uses Excel provided frameRange');
end

% parameters calculated based on user-supplied parameters
averageBacterialWidthInPixel= AVERAGEBACTERIAWIDTH/p.micronsPerPixel; % Unused
EXTRAPOLATIONLENGTH = round(1.9*averageBacterialWidthInPixel); % Maximum size of skeleton end which gets extrapolated in pixels - Independant of pixel size
paddingsize = round(averageBacterialWidthInPixel*4); % Unused

if isfield(p,'saveSummaryPlots')
    SAVESUMMARYPLOTS = p.saveSummaryPlots;
else
    if ~exist('SAVESUMMARYPLOTS','var')
        SAVESUMMARYPLOTS = 0;
    end
end

if isfield(p,'extraOutput')
    extraOutput = p.extraOutput;
else
    if ~exist('extraOutput','var')
        extraOutput = 0;
    end
end

% if isfield(p,'extraoutput')
%     % Plot with outline of all bacteria and extended skeletons
%     plot() 
%     saveas([p.analysisDir '/lengthNick/' num2str(frame) num2str(cellno) '.tif'])
% end

%% Prepare output parameters.
lastFrame = frameRange(end);
% lengths
allLengthsOfBacteriaInPixels    = cell(1,lastFrame);
allLengthsOfBacteriaInMicrons   = cell(1,lastFrame);
allLengthsOfBacteriaInPixelsMW    = cell(1,lastFrame);
allLengthsOfBacteriaInMicronsMW   = cell(1,lastFrame);
% Skeleton, area, additional data
allPixelAreaOfBacterium         = cell(1,lastFrame);    
allSkeletonXYpoleToPole         = cell(1,lastFrame);
allMinX                         = cell(1,lastFrame);
allMinY                         = cell(1,lastFrame);
allEdges                        = cell(1,lastFrame);
alldistanceAlongSkeletonPixels     = cell(1,lastFrame);
allextrapolatedDistanceEndsPixels  = cell(1,lastFrame);
allextrapolatedDistanceEndsMicrons = cell(1,lastFrame);
allExtendedSkeletons                = cell(1,lastFrame);

%% Loop over frames of this dataset
for framenr = frameRange
    
    disp(['Analyzing skeleton for frame ' num2str(framenr) ' (highest framenr =' num2str(lastFrame) ').']);

    %% % Load data for current frame of the dataset
    %e.g. load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg337.mat'
    load ([p.segmentationDir p.movieName 'seg' sprintf('%03d',framenr) '.mat']);
        % Important contents are Lc and Xreg, which respectively hold the
        % checked segmentation, and the fluorescence image, 
        
    %% Prepare loop 
    % get unique cellnos
    nonZeroIndices = (Lc(:)>0);
    allCellnos = transpose(unique(Lc(nonZeroIndices)));
    % prepare output parameters for this frame
    lengthOfBacteriaInPixelsInThisFrame = NaN(1,numel(allCellnos));
    lengthOfBacteriaInMicronsInThisFrame = NaN(1,numel(allCellnos));    
    lengthOfBacteriaInPixelsInThisFrameMW = NaN(1,numel(allCellnos));
    lengthOfBacteriaInMicronsInThisFrameMW = NaN(1,numel(allCellnos));    
    pixelAreaOfBacteriumInThisFrame =  NaN(1,numel(allCellnos));
    skeletonXYpoleToPoleInThisFrame = cell(1,numel(allCellnos));
    minXThisFrame = NaN(1,numel(allCellnos));
    minYThisFrame = NaN(1,numel(allCellnos));
    edgesThisFrame = cell(1,numel(allCellnos));
    distanceAlongSkeletonPixelsThisFrame = cell(1,numel(allCellnos));
    extrapolatedDistancePixelsEndsThisFrame = NaN(2,numel(allCellnos)); 
    extrapolatedDistanceMicronsEndsThisFrame = NaN(2,numel(allCellnos)); 
    extendedSkeletons = cell(1,numel(allCellnos));
    
    %% Loop over all cell numbers in this frame
    for cellno = allCellnos
        %% Convert one cell to x,y coordinates.        
        [y,x] = find(Lc == cellno);
        
        if extraOutput
            % show original
            h1=figure(1); clf; 
            imshow(Lc,[]);
            % show conversion
            h2=figure(2); clf;            
            plot(x,y,'.');
            axis equal;
        end
        %% % Select ROI and make image binary
        % administration required to select ROI - excludes surroundings
        minY = min(y);
        minX = min(x);
        sizeY = max(y)-minY+1;
        sizeX = max(x)-minX+1;
        % create zero array size of bacterium
        zer = zeros(sizeY,sizeX);
        % use x,y coordinates to fill it
        for framen=1:length(x)
            zer(y(framen)-minY+1,x(framen)-minX+1)=1;
        end
        binaryImage = zer;
        binaryImage = bwmorph(binaryImage,'fill'); % To avoid sporadious skeletonize error (which results in a small loop at the 'end')

        % % add padding to image (to avoid filters "seeing" edges)
        % binaryImage=padarray(binaryImage,[paddingsize,paddingsize]);

        if extraOutput
            h3=figure(3); clf;
            imshow(binaryImage);
        end
        % calculate total area of bacterium
        pixelAreaOfBacterium = sum(binaryImage(:));
        %% % Skeletonizes image
        binaryImageSkeletonized = bwmorph(binaryImage,'skel',Inf);
        %binaryImageSkeletonized = voronoiSkel(binaryImage); % downloaded this, but to tricky to get to work
        %binaryImageSkeletonized = skeleton(binaryImage); % downloaded this, but to tricky to get to work

        if extraOutput
            h4=figure(4); clf;
            imshow(binaryImageSkeletonized)
            imshow((binaryImage+binaryImageSkeletonized)/2,[])
        end
        %% % Finds edges - Gives boundary of the segmented cell
        edgesBinary = binaryImage-bwmorph(binaryImage,'erode');
        
        if extraOutput
            figure(5); clf;
            imshow(edgesBinary+binaryImageSkeletonized)
        end
        
        %% % Gets x & y values of the segmented edges, and plots them
        boundaries = bwboundaries(edgesBinary,8);
        edge = boundaries{1,1}; % Extracts correct edge of the 2 determined boundaries
        
        if extraOutput
            h71=figure(6); clf; hold on;
            axis equal;
            
            plot(edge(:,1),edge(:,2))            
            %plot(skeletonXYpoleToPole(:,1),skeletonXYpoleToPole(:,2))            
        end
        %% % Finds endings of the skeleton
        endsBeforeSpurring = bwmorph(binaryImageSkeletonized,'endpoints');
        
        if extraOutput
            h50=figure(7); clf;
            imshow(endsBeforeSpurring)
        end
        %% % Calculate number of ends before removing branches
        numEndsBeforeSpurring = sum(sum(endsBeforeSpurring,2));
        %% Finds x & y values of endings before removing branches
        k=1;        
        xyEndsBeforeSpurring = zeros(1,2,numEndsBeforeSpurring);
        for framen=1:sizeX
            for j=1:sizeY
                if endsBeforeSpurring(j,framen)==1
                    xyEndsBeforeSpurring(:,:,k)=[j,framen];
                    k=k+1;
                end
            end
        end
        
        if extraOutput
            numEndsBeforeSpurring
            xyEndsBeforeSpurring
        end
        %% % Disconnects at branch points - Not used anywhere in the script
        binaryImageDisconnected = binaryImageSkeletonized-bwmorph(binaryImageSkeletonized,'branchpoints');
        if extraOutput
            figure(8); clf;
            imshow(binaryImageDisconnected);
        end
        %% % XXX Just to try other fitting - definitely doesn't work for filamented cells XXX
        
        %{
        [xx,yy] = find(binaryImageSkeletonized==1);
        
        func = csaps(xx,yy);
        extrapolatedSpline1 = fnxtr(func);
        
        if extraOutput
            figure(71); hold on;
            axis equal;
            
            plot(edge(:,1),edge(:,2))
            
                        
            extrapolatedSpline1
            plot(xx,yy,'-')
            figure()
            fnplt(extrapolatedSpline1)            
            
        end
        %}
        
        % BWfit=fit(xx,yy,'poly9')
        % plot(BWfit)
        % BWspline=spline(xx,yy)
        %% % Removes side-branches
        % Main idea is to trim branhes until only two end-points are left
        % such that there's a branchless skeleton (i.e. the main branch).
        countSpurring=0;
        numEnds=numEndsBeforeSpurring;
        % make all pixels at the edge 0, as this otherwise can lead to
        % issues with detecting branches that are at the edge.
        binaryImageSkeletonized(1,:)   = 0;
        binaryImageSkeletonized(end,:) = 0;
        binaryImageSkeletonized(:,1)   = 0;
        binaryImageSkeletonized(:,end) = 0;
        % Continue removing spur pixels until we have branchless skeleton
        while numEnds > 2
                binaryImageSkeletonized = bwmorph(binaryImageSkeletonized,'spur');
                binaryImageSkeletonized = bwmorph(binaryImageSkeletonized,'skel'); % To prevent issue with 4-way crossings
                countSpurring = countSpurring+1;
                ends = bwmorph(binaryImageSkeletonized,'endpoints');
                numEnds = sum(sum(ends,2));
                if countSpurring>1000
                    figure(9); imshow(binaryImageSkeletonized,[]);
                    error(['Error spurring, cellnum=' num2str(cellno) ', showing current skeleton.']);                    
                end;
        end
        binarySkeletonBranchless = binaryImageSkeletonized; % Branchless skeleton
        
        if extraOutput
            countSpurring
            numEnds
        end
        %% % Finds endings of branchless skeleton
        if numEnds==1
            warning(['Skeleton is only 1 px in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
        end  
        ends = bwmorph(binarySkeletonBranchless,'endpoints');
        
        if extraOutput
            h51=figure(10); clf;
            imshow(ends)
        end
        %% % Finds x & y values of endings after removing branches
        l=1;
        xyEnds=zeros(1,2,numEnds);
        for framen=1:sizeX
            for j=1:sizeY
                if ends(j,framen)==1
                    xyEnds(:,:,l)=[j,framen];
                    l=l+1;
                end
            end
        end
        if extraOutput
            xyEnds
        end
        %% % If skeleton is 1 px --> Adds arbitrary point & finds x & y values of endings (after removing branches)
        if numEnds==1
            binarySkeletonBranchless(xyEnds(1)+1,xyEnds(2)-1) = 1; % Adds arbitrary point (8-connected)
            ends = bwmorph(binarySkeletonBranchless,'endpoints'); % Finds the 2 ends
            numEnds = numEnds+1; % Updates number of ends of branchless skeleton
            
            l=1;
            xyEnds=zeros(1,2,2);
            for framen=1:sizeX
                for j=1:sizeY
                    if ends(j,framen)==1
                    xyEnds(:,:,l)=[j,framen];
                    l=l+1;
                    end
                end
            end
        end
        if extraOutput
            xyEnds
        end
        %% % Gets x & y values of the branchless skeleton, and sets them in an array
        % Get "boundaries" of the line (result should make a loop around the line,
        % but from an arbitrary point)
        skeletonBoundary = bwboundaries(binarySkeletonBranchless,8); 
        skeletonBoundary = skeletonBoundary{1,1};
        % paste this loop 2x behind itself
        twotimesSkeletonBoundary = [skeletonBoundary; skeletonBoundary];
        leftEnd = xyEnds(:,:,1);
        % Now find one of the bacterial poles
        poleIndex = find(skeletonBoundary(:,1)==leftEnd(1) & skeletonBoundary(:,2)==leftEnd(2));
        % Now get skeletonXYpoleToPole(i,:)=v(i), with v(i,:)=[x(i),y(i)],
        % point i along the skeleton
        skeletonXYpoleToPole = twotimesSkeletonBoundary(poleIndex:poleIndex+round((size(skeletonBoundary,1)+1)/2-1),:); % -1 to correct for poleIndex
        
        if extraOutput
            h5=figure(11); clf;
            imshow((binarySkeletonBranchless+binaryImage)/2)
        end
        
        if extraOutput
            h71=figure(12); clf; hold on;
            axis equal;
            
            plot(edge(:,1),edge(:,2))            
            plot(skeletonXYpoleToPole(:,1),skeletonXYpoleToPole(:,2))            
        end
        
        %% % Sets length of dataset (coming from branchless skeleton) which gets extrapolated
        extrapolationLength = min(EXTRAPOLATIONLENGTH, length(skeletonXYpoleToPole)); % Ensures maximum window while cells are still straight on this interval
        % vq3 = interp1(array(1:20,1),array(1:20,2),'pchip');
        % bla=bspline(array(1:50,1),array(1:50,2));
        
        %% % If pieces to extrapolate contain only 1 unique x-value --> Swap x and y (transpose and switch rows/columns of variables) to ensure extrapolation works
        dyExtrapolation1 = max(skeletonXYpoleToPole(1:extrapolationLength,2))-min(skeletonXYpoleToPole(1:extrapolationLength,2));
        dxExtrapolation1 = max(skeletonXYpoleToPole(1:extrapolationLength,1))-min(skeletonXYpoleToPole(1:extrapolationLength,1));
        dydxExtrapolation1 = dyExtrapolation1/dxExtrapolation1; % Calculates how steep one end is
        dyExtrapolation2 = max(skeletonXYpoleToPole(end-(extrapolationLength-1):end,2))-min(skeletonXYpoleToPole(end-(extrapolationLength-1):end,2));
        dxExtrapolation2 = max(skeletonXYpoleToPole(end-(extrapolationLength-1):end,1))-min(skeletonXYpoleToPole(end-(extrapolationLength-1):end,1));
        dydxExtrapolation2 = dyExtrapolation2/dxExtrapolation2; % Calculates how steep other end is
        clear transposeReminder; % Makes sure this variable is cleared again every loop
        
        if numel(unique(skeletonXYpoleToPole(1:extrapolationLength,1))) < 2 || ... % Unique x-values left end
                numel(unique(skeletonXYpoleToPole(end-(extrapolationLength-1):end,1))) < 2  % Unique x-values right end
            skeletonXYpoleToPole(:,[1 2]) = skeletonXYpoleToPole(:,[2 1]); % Transposes:
            edge(:,[1 2]) = edge(:,[2 1]);
            xyEnds(:,[1 2],:) = xyEnds(:,[2 1],:);
            ends = ends';
            binarySkeletonBranchless = binarySkeletonBranchless';
            transposeReminder=1; % Added to be able to compensate for errors later
        elseif max(dydxExtrapolation1, dydxExtrapolation2) > max(1/dydxExtrapolation1, 1/dydxExtrapolation2) % Determines if transposing is beneficial
            skeletonXYpoleToPole(:,[1 2]) = skeletonXYpoleToPole(:,[2 1]); % Transposes:
            edge(:,[1 2]) = edge(:,[2 1]);
            xyEnds(:,[1 2],:) = xyEnds(:,[2 1],:);
            ends = ends';
            binarySkeletonBranchless = binarySkeletonBranchless';
            transposeReminder=2; % % Added to be able to compensate for errors later
        end
        
        %% Smooth the skeleton
        windowArray = [-SMOOTHELEMENTS:SMOOTHELEMENTS];
        %smoothSkeleton = NaN(size(currentSkeletonXYpoleToPole,1)-2*SMOOTHELEMENTS+1,2)
        smoothSkeletonXYpoleToPole = NaN(size(skeletonXYpoleToPole,1),2);
        nrIndicesInSkelet = size(skeletonXYpoleToPole,1);
        if nrIndicesInSkelet<(2*SMOOTHELEMENTS)
            warning('Extremely small skeleton, not smoothing.');
            smoothSkeletonXYpoleToPole=skeletonXYpoleToPole;
        else
            % 1st few elements
            for i = 1:SMOOTHELEMENTS
                plusminus = i-1;
                smoothSkeletonXYpoleToPole(i,1) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
                smoothSkeletonXYpoleToPole(i,2) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
            end
            % main piece of skeleton
            for i = SMOOTHELEMENTS+1:nrIndicesInSkelet-SMOOTHELEMENTS

                smoothSkeletonXYpoleToPole(i,1) = mean(skeletonXYpoleToPole(i+windowArray,1));
                smoothSkeletonXYpoleToPole(i,2) = mean(skeletonXYpoleToPole(i+windowArray,2));

            end
            % last few elements
            for i = nrIndicesInSkelet-SMOOTHELEMENTS+1:nrIndicesInSkelet
                plusminus = nrIndicesInSkelet-i;
                smoothSkeletonXYpoleToPole(i,1) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
                smoothSkeletonXYpoleToPole(i,2) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
            end
        end
        
        %% % Finds directionality of skeleton ends & creates variables to account for later on
        
        directionFindLength = min(EXTRAPOLATIONLENGTH, round(length(smoothSkeletonXYpoleToPole)/2)); % Finds appropriate length over which to decide the direction
        directionxEnd1 = smoothSkeletonXYpoleToPole(1+directionFindLength,1)-smoothSkeletonXYpoleToPole(1,1); % Positive = to the right, negative = to the left
        directionxEnd2 = smoothSkeletonXYpoleToPole(end,1)-smoothSkeletonXYpoleToPole(end-directionFindLength,1); % Positive = to the right, negative = to the left
        
        directionFactorLeft1 = 1; % Set standard values: - Only used in practice when a filamented cell has its ends (close to) perpendicular on each other
        directionFactorRight1 = 1;
        directionFactorLeft2 = 1;
        directionFactorRight2 = 1;
        if directionxEnd1 < 0 % Change values (that corresponds to the amount of extrapolation) depending upon which direction the skeleton end has:
            directionFactorLeft1  = DIRECTIONVALUECORRECT;
            directionFactorRight1 = DIRECTIONVALUEINCORRECT;
        elseif directionxEnd1 > 0
            directionFactorLeft1  = DIRECTIONVALUEINCORRECT;
            directionFactorRight1 = DIRECTIONVALUECORRECT;
        else
            warning(['First end has no directionality in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
        end
        if directionxEnd2 < 0
            directionFactorLeft2  = DIRECTIONVALUECORRECT;
            directionFactorRight2 = DIRECTIONVALUEINCORRECT;
        elseif directionxEnd2 > 0
            directionFactorLeft2  = DIRECTIONVALUEINCORRECT;
            directionFactorRight2 = DIRECTIONVALUECORRECT;
        else
            warning(['Second end has no directionality in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
        end
        
        %% % Extrapolates first end of the bacteria - fit is forced through the 'end' and extrapolates linearly outside data interval
        
        if ~exist('transposeReminder', 'var') % Make adjusting factors to compensate for relative longer extrapolation lengths at non-horizontal ends:
            adjustingFactor1 = min([1, cos(atan(dydxExtrapolation1))]);
            adjustingFactor2 = min([1, cos(atan(dydxExtrapolation2))]);
        elseif transposeReminder==1
            adjustingFactor1 = 1;
            adjustingFactor2 = 1;
        elseif transposeReminder==2
            adjustingFactor1 = min([1, cos(atan(1/dydxExtrapolation1))]);
            adjustingFactor2 = min([1, cos(atan(1/dydxExtrapolation2))]);
        end

        % actual extrapolation
        try
            %% Determine on which skeleton values left extrapolation should be based
            toFitXleft = smoothSkeletonXYpoleToPole(1:extrapolationLength,1);
            toFitYleft = smoothSkeletonXYpoleToPole(1:extrapolationLength,2);
            
            % Create extrapolation function
            %fitParamsleft = polyfit(toFitXleft,toFitYleft,1);
            %myExtrapFunctionLeft = @(x) fitParamsleft(1)*x+fitParamsleft(2);
            fitParamsleft = polyfit(toFitXleft,toFitYleft,2);
            myExtrapFunctionLeft = @(x) fitParamsleft(1)*x.^2+fitParamsleft(2)*x+fitParamsleft(3);
            
            % Create x-values for the extrapolated part
            if directionxEnd1 < 0
                % bacterial direction is towards left, extrapolate towards right
                toextrapolatexleft = [max(toFitXleft):max(toFitXleft)+extrapolationLength];
            else    
                % bacterial direction is towards right, extrapolate towards left
                toextrapolatexleft = [min(toFitXleft)-extrapolationLength:min(toFitXleft)];                
            end
            
            % create y-values based on those x-values
            extrapolatedValuesLeft = myExtrapFunctionLeft(toextrapolatexleft);
            
            %{
            % Nick's extrapolation
            func=csaps(toFitXleft,toFitYleft); % TODO MAYBE USE OTHER (POLY)FIT?
            extrapolatedSpline1 = fnxtr(func,2);
            % 'directionFactors' ensure the plotted extrapolation crosses the edge of the cell only one time (and thereby prevent calculation error)
            % 'adjustingFactor' ensures not two times the same extrapolation intersection is found
            extrapolatedSkeleton1 = fnplt(extrapolatedSpline1,[skeletonXYpoleToPole(1,1)-directionFactorLeft1*adjustingFactor1*(EXTRAPOLATIONLENGTH)... 
                skeletonXYpoleToPole(1,1)+directionFactorRight1*adjustingFactor1*(EXTRAPOLATIONLENGTH)]).';
            %}
            
            extrapolatedSkeleton1 = [toextrapolatexleft' extrapolatedValuesLeft'];
        catch
            cellno
            figure(13); imshow(binaryImage+binaryImageSkeletonized,[]);
            smoothSkeletonXYpoleToPole
            error('Extrapolation failed.');
        end 
        
        if extraOutput
            %%           
            
            extrapolatedSkeleton1
            %figure()            
            %fnplt(extrapolatedSpline1,[smoothSkeletonXYpoleToPole(1,1)-directionFactorLeft1*adjustingFactor1*(EXTRAPOLATIONLENGTH)... 
            %    smoothSkeletonXYpoleToPole(1,1)+directionFactorRight1*adjustingFactor1*(EXTRAPOLATIONLENGTH)])
        end
        %% % Extrapolates second end of the bacteria - fit is forced through the 'end' and extrapolates linearly outside data interval
        %% Determine on which skeleton values left extrapolation should be based
        toFitXright = smoothSkeletonXYpoleToPole(length(smoothSkeletonXYpoleToPole)-(extrapolationLength-1):length(smoothSkeletonXYpoleToPole),1);
        toFitYright = smoothSkeletonXYpoleToPole(length(smoothSkeletonXYpoleToPole)-(extrapolationLength-1):length(smoothSkeletonXYpoleToPole),2);

        % Create extrapolation function
        %fitParamsRight = polyfit(toFitXright,toFitYright,1);
        %myExtrapFunctionRight = @(x) fitParamsRight(1)*x+fitParamsRight(2);
        fitParamsRight = polyfit(toFitXright,toFitYright,2);
        myExtrapFunctionRight = @(x) fitParamsRight(1)*x.^2+fitParamsRight(2)*x+fitParamsRight(3);
        

        % Create x-values for the extrapolated part
        if directionxEnd2 < 0
            toextrapolatexright = [min(toFitXright)-extrapolationLength:min(toFitXright)];
        else
            toextrapolatexright = [max(toFitXright):max(toFitXright)+extrapolationLength];            
        end

        % create y-values based on those x-values
        extrapolatedValuesRight = myExtrapFunctionRight(toextrapolatexright);                
        
        % Nick's extrapolation (old)
        %{
        func2=csaps(skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extrapolationLength-1):length(skeletonXYpoleToPole),1),skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extrapolationLength-1):length(skeletonXYpoleToPole),2));
        extrapolatedSpline2 = fnxtr(func2);
        
        
        % 'directionFactors' ensure the plotted extrapolation crosses the edge of the cell only one time (and thereby prevent calculation error)
        % 'adjustingFactor' ensures not two times the same extrapolation intersection is found
        extrapolatedSkeleton2 = fnplt(extrapolatedSpline2,[skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)-directionFactorLeft2*adjustingFactor2*(EXTRAPOLATIONLENGTH)... 
            skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)+directionFactorRight2*adjustingFactor2*(EXTRAPOLATIONLENGTH)]).';
        %}
        
        extrapolatedSkeleton2 = [toextrapolatexright' extrapolatedValuesRight'];
        
        if extraOutput
            extrapolatedSkeleton2
            %figure()            
            %fnplt(extrapolatedSpline2,[smoothSkeletonXYpoleToPole(length(smoothSkeletonXYpoleToPole),1)-directionFactorLeft2*adjustingFactor2*(EXTRAPOLATIONLENGTH)... 
            %    smoothSkeletonXYpoleToPole(length(smoothSkeletonXYpoleToPole),1)+directionFactorRight2*adjustingFactor2*(EXTRAPOLATIONLENGTH)])
        end
        %% % Plot extrapolations and segmented edges
        if extraOutput
            %%
            h72=figure(14); clf; axis equal; hold on;            
            
            plot(edge(:,1),edge(:,2));
            %plot(extrapolatedSkeleton1(:,1),extrapolatedSkeleton1(:,2)); % Nick's extrapolation (old)
            %plot(extrapolatedSkeleton2(:,1),extrapolatedSkeleton2(:,2)); % Nick's extrapolation (old)
            plot(smoothSkeletonXYpoleToPole(:,1),smoothSkeletonXYpoleToPole(:,2));
            
            % right end
            plot(toFitXleft,toFitYleft,'.','MarkerSize',10);
            plot(toextrapolatexleft,extrapolatedValuesLeft,'-k');
            % left end
            plot(toFitXright,toFitYright,'.','MarkerSize',10);
            plot(toextrapolatexright,extrapolatedValuesRight,'-k');
        end
        %% % Determine intersection point and with that the correction length for one end
        % Create parameter arrays
        disx=zeros(length(extrapolatedSkeleton1),length(edge));
        disy=zeros(length(extrapolatedSkeleton1),length(edge));
        distot=zeros(length(extrapolatedSkeleton1),length(edge));
        
        for framen=1:length(extrapolatedSkeleton1)
            for j=1:length(edge)
                disx(framen,j) = edge(j,1)-extrapolatedSkeleton1(framen,1); % Distance in x between every point on the edge and the extrapolation
                disy(framen,j) = edge(j,2)-extrapolatedSkeleton1(framen,2); % Distance in y between every point on the edge and the extrapolation
                %distot(framen,j) = disx(framen,j)+disy(framen,j); % Distance cityblock-way
                distot(framen,j) = sqrt(disx(framen,j).^2+disy(framen,j).^2); % Distance Pythagoras-way instead of cityblock-way
            end
        end        
        minimumDistance1 = min(min(distot)); % Finds minimum of distance matrix --> Shortest distance between edge and extrapolation
        if minimumDistance1 > EXTRAPOLATIONLENGTH*ERRORINTERSECT % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on first end might have gone wrong in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            minimumDistance1
        elseif minimumDistance1 > EXTRAPOLATIONLENGTH*ERRORINTERSECTWARNING % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on first end is just above error threshold in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            minimumDistance1
        end 
        % Find coördinates corresponding to closest points between edge and extrapolation
        for framen=1:length(extrapolatedSkeleton1)
            for j=1:length(edge)
                if distot(framen,j)==minimumDistance1
                    icor=framen; % Saves which specific point of the extrapolation was closest
                    jcor=j; % Saves which specific point of the edge was closest
                end
            end
        end
        
        %extrapolatedIntersection1 = extrapolatedSkeleton1(icor,:);
        extrapolatedIntersection1 = edge(jcor,:); % Point on the cell boundary that gets extrapolated to        
        
        % create the x,y coordinates for the extrapolated piece within cell
        % boundaries
        if directionxEnd1<0
            extrapolatedSkeleton1WithinEdge = extrapolatedSkeleton1(icor:-1:1,:);
        else
            extrapolatedSkeleton1WithinEdge = extrapolatedSkeleton1(icor:end,:);            
        end        
        
        % determine distance to both ends from found intersection point
        extrapolatedDistance1FirstEnd = pdist2(extrapolatedIntersection1,xyEnds(1,:,1));
        extrapolatedDistance1SecondEnd = pdist2(extrapolatedIntersection1,xyEnds(1,:,numEnds));
        % take smallest as relevant extrapolation distance of the first end        
        extrapolatedDistance1 = min([extrapolatedDistance1FirstEnd extrapolatedDistance1SecondEnd]);        
        if extrapolatedDistance1 > ERRORFACTOREXTRAPOLATIONLENGTH*EXTRAPOLATIONLENGTH % Warns when the extrapolated distance becomes quite large
            warning(['Extrapolation is quite big on the first end in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            extrapolatedDistance1
        end
        
        if extraOutput
            extrapolatedIntersection1
            minimumDistance1
            edge(jcor,:)            
            xyEnds(1,:,1)
            extrapolatedDistance1            
        end
        
        %% % Determine other intersection point and with that the correction length for the other end
        % Create parameter arrays
        disx2=zeros(length(extrapolatedSkeleton2),length(edge));
        disy2=zeros(length(extrapolatedSkeleton2),length(edge));
        distot2=zeros(length(extrapolatedSkeleton2),length(edge));
        
        for framen=1:length(extrapolatedSkeleton2)
            for j=1:length(edge)
                disx2(framen,j) = edge(j,1)-extrapolatedSkeleton2(framen,1); % Distance in x between every point on the edge and the extrapolation
                disy2(framen,j) = edge(j,2)-extrapolatedSkeleton2(framen,2); % Distance in y between every point on the edge and the extrapolation
                %distot2(framen,j) = disx2(framen,j)+disy2(framen,j); % Distance cityblock-way
                distot2(framen,j) = sqrt(disx2(framen,j).^2+disy2(framen,j).^2); % Distance Pythagoras-way instead of cityblock-way
            end
        end        
        minimumDistance2 = min(min(distot2)); % Finds minimum of second distance matrix --> Shortest distance between edge and extrapolation
        if minimumDistance2 > EXTRAPOLATIONLENGTH*ERRORINTERSECT % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on second end might have gone wrong in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            minimumDistance2
        elseif minimumDistance2 > EXTRAPOLATIONLENGTH*ERRORINTERSECTWARNING % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on second end is just above error threshold in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            minimumDistance2
        end
        % Find coördinates corresponding to closest points between edge and extrapolation
        for framen=1:length(extrapolatedSkeleton2)
            for j=1:length(edge)
                if distot2(framen,j)==minimumDistance2
                    icor2=framen; % Saves which specific point of the extrapolation was closest
                    jcor2=j; % Saves which specific point of the edge was closest
                end
            end
        end
           
        %extrapolatedIntersection2 = extrapolatedSkeleton2(icor2,:);
        extrapolatedIntersection2 = edge(jcor2,:); % Second point on the cell boundary that gets extrapolated to
        
        % create the x,y coordinates for the extrapolated piece within cell
        % boundaries
        if directionxEnd2<0
            extrapolatedSkeleton2WithinEdge = extrapolatedSkeleton2(end:-1:icor2,:);
        else            
            extrapolatedSkeleton2WithinEdge = extrapolatedSkeleton2(1:icor2,:);
        end   
        
        % determine distance to both ends from found intersection point
        extrapolatedDistance2FirstEnd = pdist2(extrapolatedIntersection2,xyEnds(1,:,1));
        extrapolatedDistance2SecondEnd = pdist2(extrapolatedIntersection2,xyEnds(1,:,numEnds));
        % take smallest as relevant extrapolation distance of the second end
        extrapolatedDistance2 = min([extrapolatedDistance2FirstEnd extrapolatedDistance2SecondEnd]);
        if extrapolatedDistance2 > ERRORFACTOREXTRAPOLATIONLENGTH*EXTRAPOLATIONLENGTH % Warns when the extrapolated distance becomes quite large
            warning(['Extrapolation is quite big on the second end in frame ' num2str(framenr) ' cell ' num2str(cellno)]);
            extrapolatedDistance2
        end
                
        if extraOutput
            extrapolatedIntersection2
            minimumDistance2
            edge(jcor2,:)            
            xyEnds(1,:,numEnds)
            extrapolatedDistance2
        end
        
        
        %% % Additional check for errors
        if isequal(extrapolatedIntersection1, extrapolatedIntersection2) && length(smoothSkeletonXYpoleToPole) > ERRORVALUESKELETONLENGTH
            warning(['Extrapolation distances are the same in frame ' num2str(framenr) ' cell ' num2str(cellno) ' with skeletonlength ' num2str(length(smoothSkeletonXYpoleToPole))]);
            disp(extrapolatedDistance1);
        end
        %% % Calculates length of branchless skeleton, and the total estimated length (by adding the extrapolated lengths of the ends) 
        distanceMask = ends; % Binary array with the ends of the branchless skeleton indicated as 1's
        extractEnd = xyEnds(1,:,1); % Get x and y value of the first end
        distanceMask(extractEnd(1),extractEnd(2)) = 0; % Set this found end to a 0 in the array
        D = bwdistgeodesic(binarySkeletonBranchless,distanceMask,'quasi-euclidean'); % Computes distance along the branchless skeleton

        % Writes the computed distances (along the branchless skeleton) in a smaller array
        distanceAlongSkeletonPixels = D(sub2ind(size(D),round(smoothSkeletonXYpoleToPole(:,1)),round(smoothSkeletonXYpoleToPole(:,2))));
        
        % Get end-to-end distance along the branchless skeleton
        EndToEndDistanceSkeleton = max(max(D)); % Same as "max(distanceAlongSkeletonPixels)"

        if extraOutput
            EndToEndDistanceSkeleton
        end 
        
        % now before exporting backtranspose if necessary
        if exist('transposeReminder','var') 
            smoothSkeletonXYpoleToPole(:,[1 2]) = smoothSkeletonXYpoleToPole(:,[2 1]); % Transposes:
            edge(:,[1 2]) = edge(:,[2 1]);
            xyEnds(:,[1 2],:) = xyEnds(:,[2 1],:);
            ends = ends';
            binarySkeletonBranchless = binarySkeletonBranchless';
            
            extrapolatedSkeleton1(:,[1 2]) = extrapolatedSkeleton1(:,[2 1]);
            extrapolatedSkeleton2(:,[1 2]) = extrapolatedSkeleton2(:,[2 1]);
            
            extrapolatedSkeleton1WithinEdge(:,[1 2]) = extrapolatedSkeleton1WithinEdge(:,[2 1]);
            extrapolatedSkeleton2WithinEdge(:,[1 2]) = extrapolatedSkeleton2WithinEdge(:,[2 1]);

        end
        
        %% Length of skeleton pieces
        mainPartLength = sum(MW_distancesbetweenlinesegments(smoothSkeletonXYpoleToPole(:,1)',smoothSkeletonXYpoleToPole(:,2)'));
        extrapolatedPart1Length = sum(MW_distancesbetweenlinesegments(extrapolatedSkeleton1WithinEdge(:,1)',extrapolatedSkeleton1WithinEdge(:,2)'));
        extrapolatedPart2Length = sum(MW_distancesbetweenlinesegments(extrapolatedSkeleton2WithinEdge(:,1)',extrapolatedSkeleton2WithinEdge(:,2)'));
        totalMWLengthPixels = mainPartLength + extrapolatedPart1Length + extrapolatedPart2Length;
        
         %% summary plot
        if extraOutput || SAVESUMMARYPLOTS
            h72 = figure(72); clf; axis equal; hold on;
            if SAVESUMMARYPLOTS % hide figure if it is saved to disk
                set(h72,'visible','off');
            end
            
            plot(edge(:,1),edge(:,2));
            %plot(extrapolatedSkeleton1(:,1),extrapolatedSkeleton1(:,2)); % Nick's extrapolation (old)
            %plot(extrapolatedSkeleton2(:,1),extrapolatedSkeleton2(:,2)); % Nick's extrapolation (old)
            plot(smoothSkeletonXYpoleToPole(:,1),smoothSkeletonXYpoleToPole(:,2));
            
            % left end
            %plot(toFitXleft,toFitYleft,'.','MarkerSize',10);
            plot(extrapolatedSkeleton1WithinEdge(:,1),extrapolatedSkeleton1WithinEdge(:,2),'-ko');
            % right end
            %plot(toFitXright,toFitYright,'.','MarkerSize',10);
            plot(extrapolatedSkeleton2WithinEdge(:,1),extrapolatedSkeleton2WithinEdge(:,2),'-ko');
            
            if exist('transposeReminder','var') 
                plot(toFitYleft,toFitXleft,'.','MarkerSize',10);
                plot(toFitYright,toFitXright,'.','MarkerSize',10);
            else
                plot(toFitXleft,toFitYleft,'.','MarkerSize',10);
                plot(toFitXright,toFitYright,'.','MarkerSize',10);
            end
            
            % Save the figure if desired
            if SAVESUMMARYPLOTS
                theOutputdir = [p.analysisDir 'straightenedCells\skeletons\'];

                if ~exist(theOutputdir,'dir')
                    mkdir(theOutputdir)
                end;

                saveas(h72,[theOutputdir 'fr' num2str(framenr) 'cellno' num2str(cellno) '.tif']);
            end
        end
        
        %% export length data 
        lengthOfBacteriaInPixelsInThisFrame(cellno)  = EndToEndDistanceSkeleton+extrapolatedDistance1+extrapolatedDistance2;
        lengthOfBacteriaInMicronsInThisFrame(cellno) = lengthOfBacteriaInPixelsInThisFrame(cellno)*p.micronsPerPixel;
        
        lengthOfBacteriaInPixelsInThisFrameMW(cellno) = totalMWLengthPixels;
        lengthOfBacteriaInMicronsInThisFrameMW(cellno) = totalMWLengthPixels*p.micronsPerPixel;
        % export additional data
        pixelAreaOfBacteriumInThisFrame(cellno) = pixelAreaOfBacterium;
        skeletonXYpoleToPoleInThisFrame{cellno} = skeletonXYpoleToPole;
        smoothSkeletonXYpoleToPoleInThisFrame{cellno} = smoothSkeletonXYpoleToPole;
        minXThisFrame(cellno) = minX;
        minYThisFrame(cellno) = minY;
        edgesThisFrame{cellno} = edge;
        distanceAlongSkeletonPixelsThisFrame{cellno} = distanceAlongSkeletonPixels;
        extrapolatedDistancePixelsEndsThisFrame(:,cellno) = [extrapolatedDistance1 extrapolatedDistance2];
        extrapolatedDistanceMicronsEndsThisFrame(:,cellno) = [extrapolatedDistance1 extrapolatedDistance2]*p.micronsPerPixel;
        % Skeleton XY values + extrapolated parts
        extendedSkeletons{cellno} = [extrapolatedSkeleton1WithinEdge; smoothSkeletonXYpoleToPole; extrapolatedSkeleton2WithinEdge];
        
    end
    
    %% Creates summary information
    % lengths
    allLengthsOfBacteriaInPixels{framenr} = lengthOfBacteriaInPixelsInThisFrame;
    allLengthsOfBacteriaInMicrons{framenr} = lengthOfBacteriaInMicronsInThisFrame;
    
    allLengthsOfBacteriaInPixelsMW{framenr} = lengthOfBacteriaInPixelsInThisFrameMW;
    allLengthsOfBacteriaInMicronsMW{framenr} =lengthOfBacteriaInMicronsInThisFrameMW;
    
    % Skeleton, area, additional data
    allPixelAreaOfBacterium{framenr} = pixelAreaOfBacteriumInThisFrame;    
    allsmoothSkeletonXYpoleToPole{framenr} = smoothSkeletonXYpoleToPoleInThisFrame;
    allSkeletonXYpoleToPole{framenr} = skeletonXYpoleToPoleInThisFrame;    
    allMinX{framenr} = minXThisFrame;
    allMinY{framenr} = minYThisFrame;
    allEdges{framenr} = edgesThisFrame;
    alldistanceAlongSkeletonPixels{framenr} = distanceAlongSkeletonPixelsThisFrame;
    allextrapolatedDistanceEndsPixels{framenr} = extrapolatedDistancePixelsEndsThisFrame;
    allextrapolatedDistanceEndsMicrons{framenr} = extrapolatedDistanceMicronsEndsThisFrame;
    
    allExtendedSkeletons{framenr} = extendedSkeletons;
    
    %% Saves summary information
    save([p.tracksDir p.movieName '-skeletonData.mat'],...
        'allLengthsOfBacteriaInPixels','allLengthsOfBacteriaInMicrons',...
        'allLengthsOfBacteriaInPixelsMW','allLengthsOfBacteriaInMicronsMW',...
        'allPixelAreaOfBacterium','allSkeletonXYpoleToPole',...
        'allsmoothSkeletonXYpoleToPole',...
        'allMinX','allMinY',...
        'allEdges','alldistanceAlongSkeletonPixels',...
        'allextrapolatedDistanceEndsPixels','allextrapolatedDistanceEndsMicrons','allExtendedSkeletons');    
end

%{ 
%% Old code

% int2=zeros(1,length(points));
% for i=1:length(points)
%     int2(1,i)=min(disx(i,:))+min(disy(i,:));
% end
% int2;
% %
% Alternative script to get rid of side-branches
% skel= bwmorph(binaryImage,'skel',Inf);
% skel = testBacterium>0;
% 
% B = bwmorph(skel, 'branchpoints');
% E = bwmorph(skel, 'endpoints');
% 
% [y,x] = find(E);
% B_loc = find(B);
% 
% Dmask = false(size(skel));
% for k = 1:numel(x)
%     D = bwdistgeodesic(skel,x(k),y(k));
%     distanceToBranchPt = min(D(B_loc));
%     Dmask(D < distanceToBranchPt) =true;
% end
% skelD = skel - Dmask;
% h60=figure(); clf;
% imshow(skelD);
% hold all;
% [y,x] = find(B); plot(x,y,'ro')
% leng=sum(sum(skel,2))

%% Plot the total skeleton on top of the bacterium
%}

%{
% %% smooth using circular neighborhood
        % % Seems to make it a little better for my cells
        % SCALING = 1; % scaling is not so useful actually, 1 means no scaling is applied
        % 
        % work_img = binaryImage;
        % 
        % % rescaling, maybe higher resolution better?
        % work_img = imresize(work_img,SCALING,'nearest');
        % 
        % % set circular neighborhood
        % % note that the size of the disk is pretty critical, and due to the size 
        % % of the bacterium is too small to have an actual effect
        % sizeOfDisk=round(averageBacterialWidthInPixel)*.5*SCALING;
        % se=strel('disk',sizeOfDisk);
        % 
        % % apply filters
        % %work_img = imerode(work_img,se);
        % %work_img = imdilate(work_img,se);
        % work_img = imclose(work_img,se); % works best
        % 
        % % show results
        % figure(); 
        % subplot(1,2,1);
        % imshow(binaryImage);
        % subplot(1,2,2)
        % imshow(work_img);
        % 
        % if 0 % to turn on/off this filter
        %     binaryImage = work_img;
        % end
%}