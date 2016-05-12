function [p, allLengthsOfBacteriaInPixels, allLengthsOfBacteriaInMicrons] = NDL_lengthforfillamentedcells(p, frameRange) 
%function [p,schnitzcells] = NDL_lengthforfillamentedcells(p) 
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
% p.extraOutput    if this parameter is set, extra output will be shown in
%                   plots and in command window.
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
% >> NDL_lengthforfillamentedcells(p, settings.frameRangeFull) 
%
%

% function parameters set by user
AVERAGEBACTERIAWIDTH = .5; % In micron
% EXTRAPOLATIONLENGTH = 30; % In pixels - Dependant of pixel size --> CHANGE WHEN PIXEL SIZE IS DIFFERENT

% frameRange = unique([schnitzcells(:).frame_nrs]);
if exist('settings', 'var') == 1
    frameRange = settings.frameRangeFull; % Sets framerange to the full framerange provided in the Excel file
    warning('Uses Excel provided frameRange');
end

% parameters calculated based on user-supplied parameters
averageBacterialWidthInPixel= AVERAGEBACTERIAWIDTH/p.micronsPerPixel; % Unused
EXTRAPOLATIONLENGTH = round(1.8*averageBacterialWidthInPixel)+1; % Maximum size of skeleton end which gets extrapolated in pixels - Independant of pixel size
paddingsize = round(averageBacterialWidthInPixel*4); % Unused

if isfield(p,'extraOutput')
    extraOutput = p.extraOutput;
else
    extraOutput = 0;
end

% if isfield(p,'extraoutput')
%     % Plot with outline of all bacteria and extended skeletons
%     plot() 
%     saveas([p.analysisDir '/lengthNick/' num2str(frame) num2str(cellno) '.tif'])
% end
%% Loop over frames of this dataset
% Prepare output parameters.
lastFrame = frameRange(end);
% lengths
allLengthsOfBacteriaInPixels    = cell(1,lastFrame);
allLengthsOfBacteriaInMicrons   = cell(1,lastFrame);
% Skeleton, area, additional data
allPixelAreaOfBacterium         = cell(1,lastFrame);    
allSkeletonXYpoleToPole         = cell(1,lastFrame);
allMinX                         = cell(1,lastFrame);
allMinY                         = cell(1,lastFrame);
allEdges                        = cell(1,lastFrame);
alldistanceAlongSkeletonPixels     = cell(1,lastFrame);
allextrapolatedDistanceEndsPixels  = cell(1,lastFrame);
allextrapolatedDistanceEndsMicrons = cell(1,lastFrame);


for framenr = frameRange
    
    disp(['Analyzing frame ' num2str(framenr) ' (highest framenr =' num2str(lastFrame) ').']);

    %% % Load data for current frame of the dataset
    %e.g. load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg337.mat'
    load ([p.segmentationDir p.movieName 'seg' sprintf('%03d',framenr) '.mat']);
        % Important contents are Lc and Xreg, which respectively hold the
        % checked segmentation, and the fluorescence image, 
    %% Loop over all cell numbers in this frame
    % get unique cellnos
    nonZeroIndices = (Lc(:)>0);
    allCellnos = transpose(unique(Lc(nonZeroIndices)));
    % prepare output parameters for this frame
    lengthOfBacteriaInPixelsInThisFrame = NaN(1,numel(allCellnos));
    lengthOfBacteriaInMicronsInThisFrame = NaN(1,numel(allCellnos));    
    pixelAreaOfBacteriumInThisFrame =  NaN(1,numel(allCellnos));
    skeletonXYpoleToPoleInThisFrame = cell(1,numel(allCellnos));
    minXThisFrame = NaN(1,numel(allCellnos));
    minYThisFrame = NaN(1,numel(allCellnos));
    edgesThisFrame = cell(1,numel(allCellnos));
    distanceAlongSkeletonPixelsThisFrame = cell(1,numel(allCellnos));
    extrapolatedDistancePixelsEndsThisFrame = NaN(2,numel(allCellnos)); 
    extrapolatedDistanceMicronsEndsThisFrame = NaN(2,numel(allCellnos)); 
    % loop
    for cellnum = allCellnos
        %% Convert one cell to x,y coordinates.        
        [y,x] = find(Lc == cellnum);
        
        if extraOutput
            % show original
            figure(1); clf; 
            imshow(Lc,[]);
            % show conversion
            figure(2); clf;
            axis equal;
            plot(x,y,'.');
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
            figure(3);
            imshow(binaryImage);
        end
        % calculate total area of bacterium
        pixelAreaOfBacterium = sum(binaryImage(:));
        %% % Skeletonizes image
        binaryImageSkeletonized = bwmorph(binaryImage,'skel',Inf);
        %binaryImageSkeletonized = voronoiSkel(binaryImage); % downloaded this, but to tricky to get to work
        %binaryImageSkeletonized = skeleton(binaryImage); % downloaded this, but to tricky to get to work

        if extraOutput
            figure(4); clf;
            imshow(binaryImageSkeletonized)
            imshow((binaryImage+binaryImageSkeletonized)/2,[])
        end
        %% % Finds edges - Gives boundary of the segmented cell
        edgesBinary = binaryImage-bwmorph(binaryImage,'erode');
        
        if extraOutput
            figure(); clf;
            imshow(edgesBinary+binaryImageSkeletonized)
        end
        %% % Finds endings of the skeleton
        endsBeforeSpurring = bwmorph(binaryImageSkeletonized,'endpoints');
        
        if extraOutput
            figure(50); clf;
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
            figure(); clf;
            imshow(binaryImageDisconnected);
        end
        %% % XXX Just to try other fitting - definitely doesn't work for filamented cells XXX
        [xx,yy] = find(binaryImageSkeletonized==1);
        
        func = csaps(xx,yy);
        extrapolatedSpline1 = fnxtr(func);
        
        if extraOutput
            extrapolatedSpline1
            plot(xx,yy,'.')
            figure()
            fnplt(extrapolatedSpline1)
        end
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
                    figure(); imshow(binaryImageSkeletonized,[]);
                    error(['Error spurring, cellnum=' num2str(cellnum) ', showing current skeleton.']);                    
                end;
        end
        binarySkeletonBranchless = binaryImageSkeletonized; % Branchless skeleton
        
        if extraOutput
            countSpurring
            numEnds
        end
        %% % Finds endings of branchless skeleton
        if numEnds==1
            warning(['Skeleton is only 1 px in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
        end  
        ends = bwmorph(binarySkeletonBranchless,'endpoints');
        
        if extraOutput
            figure(51); clf;
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
            figure(5); clf;
            imshow((binarySkeletonBranchless+binaryImage)/2)
        end
        %% % Gets x & y values of the segmented edges, and plots them
        boundaries = bwboundaries(edgesBinary,8);
        edge = boundaries{1,1}; % Extracts correct edge of the 2 determined boundaries
        
        if extraOutput
            figure(71)
            plot(edge(:,1),edge(:,2))
            hold on
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
        %% % Extrapolates first end of the bacteria - fit is forced through the 'end' and extrapolates linearly outside data interval
        if ~exist('transposeReminder', 'var') % Make adjusting factors to compensate for relative longer extrapolation lengths at non-horizontal ends
            adjustingFactor1 = min([1, cos(atan(dydxExtrapolation1))]);
            adjustingFactor2 = min([1, cos(atan(dydxExtrapolation2))]);
        elseif transposeReminder==1
            adjustingFactor1 = 1;
            adjustingFactor2 = 1;
        elseif transposeReminder==2
            adjustingFactor1 = min([1, cos(atan(1/dydxExtrapolation1))]);
            adjustingFactor2 = min([1, cos(atan(1/dydxExtrapolation2))]);
        end

        try
            func=csaps(skeletonXYpoleToPole(1:extrapolationLength,1),skeletonXYpoleToPole(1:extrapolationLength,2)); % TODO MAYBE USE OTHER (POLY)FIT?
            extrapolatedSpline1 = fnxtr(func,2);
            % 'countSpurring' ensures the plotted extrapolation crosses the edge of the cell for cells with small branchless skeletons (many iterations of 'spur' & 'skel')
            % 'adjustingFactor' ensures not two times the same extrapolation intersection is found
            extrapolatedSkeleton1 = fnplt(extrapolatedSpline1,[skeletonXYpoleToPole(1,1)-adjustingFactor1*(EXTRAPOLATIONLENGTH)... 
                skeletonXYpoleToPole(1,1)+adjustingFactor1*(EXTRAPOLATIONLENGTH)]).';
        catch
            cellnum
            figure(); imshow(binaryImage+binaryImageSkeletonized,[]);
            skeletonXYpoleToPole
            error('Extrapolation failed.');
        end 
        
        if extraOutput
            extrapolatedSkeleton1
            figure()            
            fnplt(extrapolatedSpline1,[skeletonXYpoleToPole(1,1)-adjustingFactor1*(EXTRAPOLATIONLENGTH) skeletonXYpoleToPole(1,1)+adjustingFactor1*(EXTRAPOLATIONLENGTH)])
        end
        %% % Extrapolates second end of the bacteria - fit is forced through the 'end' and extrapolates linearly outside data interval
        func2=csaps(skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extrapolationLength-1):length(skeletonXYpoleToPole),1),skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extrapolationLength-1):length(skeletonXYpoleToPole),2));
        extrapolatedSpline2 = fnxtr(func2);
        % 'countSpurring' ensures the plotted extrapolation crosses the edge of the cell for cells with small branchless skeletons (many iterations of 'spur' & 'skel')
        % 'adjustingFactor' ensures not two times the same extrapolation intersection is found
        extrapolatedSkeleton2 = fnplt(extrapolatedSpline2,[skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)-adjustingFactor2*(EXTRAPOLATIONLENGTH)... 
            skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)+adjustingFactor2*(EXTRAPOLATIONLENGTH)]).';
        
        if extraOutput
            extrapolatedSkeleton2
            figure()            
            fnplt(extrapolatedSpline2,[skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)-adjustingFactor2*(EXTRAPOLATIONLENGTH) skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)+adjustingFactor2*(EXTRAPOLATIONLENGTH)])
        end
        %% % Plot extrapolations and segmented edges
        if extraOutput
            figure(72)
            plot(edge(:,1),edge(:,2))
            hold on
            plot(extrapolatedSkeleton1(:,1),extrapolatedSkeleton1(:,2))
            hold on
            plot(extrapolatedSkeleton2(:,1),extrapolatedSkeleton2(:,2))
            hold on
            plot(skeletonXYpoleToPole(:,1),skeletonXYpoleToPole(:,2))
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
        if minimumDistance1 > EXTRAPOLATIONLENGTH/5 % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on first end might have gone wrong in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
            minimumDistance1
        elseif minimumDistance1 > EXTRAPOLATIONLENGTH/10 % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on first end is just above error threshold in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
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
        
        % determine distance to both ends from found intersection point
        extrapolatedDistance1FirstEnd = pdist2(extrapolatedIntersection1,xyEnds(1,:,1));
        extrapolatedDistance1SecondEnd = pdist2(extrapolatedIntersection1,xyEnds(1,:,numEnds));
        % take smallest as relevant extrapolation distance of the first end
        extrapolatedDistance1 = min([extrapolatedDistance1FirstEnd extrapolatedDistance1SecondEnd]);
        if extrapolatedDistance1 > EXTRAPOLATIONLENGTH % Warns when the extrapolated distance becomes quite large
            warning(['Extrapolation is quite big on the first end in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
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
        if minimumDistance2 > EXTRAPOLATIONLENGTH/5 % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on second end might have gone wrong in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
            minimumDistance2
        elseif minimumDistance2 > EXTRAPOLATIONLENGTH/10 % Warns when found extrapolated intersection is not close to the segmented boundary
            warning(['Extrapolation on second end is just above error threshold in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
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
        
        % determine distance to both ends from found intersection point
        extrapolatedDistance2FirstEnd = pdist2(extrapolatedIntersection2,xyEnds(1,:,1));
        extrapolatedDistance2SecondEnd = pdist2(extrapolatedIntersection2,xyEnds(1,:,numEnds));
        % take smallest as relevant extrapolation distance of the second end
        extrapolatedDistance2 = min([extrapolatedDistance2FirstEnd extrapolatedDistance2SecondEnd]);
        if extrapolatedDistance2 > EXTRAPOLATIONLENGTH % Warns when the extrapolated distance becomes quite large
            warning(['Extrapolation is quite big on the second end in frame ' num2str(framenr) ' cell ' num2str(cellnum)]);
            extrapolatedDistance2
        end
                
        if extraOutput
            extrapolatedIntersection2
            minimumDistance2
            edge(jcor2,:)            
            xyEnds(1,:,numEnds)
            extrapolatedDistance2
        end
        
        %% % Calculates length of branchless skeleton, and the total estimated length (by adding the extrapolated lengths of the ends) 
        distanceMask = ends; % Binary array with the ends of the branchless skeleton indicated as 1's
        extractEnd = xyEnds(1,:,1); % Get x and y value of the first end
        distanceMask(extractEnd(1),extractEnd(2)) = 0; % Set this found end to a 0 in the array
        D = bwdistgeodesic(binarySkeletonBranchless,distanceMask,'quasi-euclidean'); % Computes distance along the branchless skeleton

        % Writes the computed distances (along the branchless skeleton) in a smaller array
        distanceAlongSkeletonPixels = D(sub2ind(size(D),round(skeletonXYpoleToPole(:,1)),round(skeletonXYpoleToPole(:,2))));
        
        % Get end-to-end distance along the branchless skeleton
        EndToEndDistanceSkeleton = max(max(D)); % Same as "max(distanceAlongSkeletonPixels)"

        if extraOutput
            EndToEndDistanceSkeleton
        end 
        % export length data 
        lengthOfBacteriaInPixelsInThisFrame(cellnum)  = EndToEndDistanceSkeleton+extrapolatedDistance1+extrapolatedDistance2;
        lengthOfBacteriaInMicronsInThisFrame(cellnum) = lengthOfBacteriaInPixelsInThisFrame(cellnum)*p.micronsPerPixel;
        % export additional data
        pixelAreaOfBacteriumInThisFrame(cellnum) = pixelAreaOfBacterium;
        skeletonXYpoleToPoleInThisFrame{cellnum} = skeletonXYpoleToPole;
        minXThisFrame(cellnum) = minX;
        minYThisFrame(cellnum) = minY;
        edgesThisFrame{cellnum} = edge;
        distanceAlongSkeletonPixelsThisFrame{cellnum} = distanceAlongSkeletonPixels;
        extrapolatedDistancePixelsEndsThisFrame(:,cellnum) = [extrapolatedDistance1 extrapolatedDistance2];
        extrapolatedDistanceMicronsEndsThisFrame(:,cellnum) = [extrapolatedDistance1 extrapolatedDistance2]*p.micronsPerPixel;
    end
    % Saves important information:    
    % lengths
    allLengthsOfBacteriaInPixels{framenr} = lengthOfBacteriaInPixelsInThisFrame;
    allLengthsOfBacteriaInMicrons{framenr} = lengthOfBacteriaInMicronsInThisFrame;
    
    % Skeleton, area, additional data
    allPixelAreaOfBacterium{framenr} = pixelAreaOfBacteriumInThisFrame;    
    allSkeletonXYpoleToPole{framenr} = skeletonXYpoleToPoleInThisFrame;
    allMinX{framenr} = minXThisFrame;
    allMinY{framenr} = minYThisFrame;
    allEdges{framenr} = edgesThisFrame;
    alldistanceAlongSkeletonPixels{framenr} = distanceAlongSkeletonPixelsThisFrame;
    allextrapolatedDistanceEndsPixels{framenr} = extrapolatedDistancePixelsEndsThisFrame;
    allextrapolatedDistanceEndsMicrons{framenr} = extrapolatedDistanceMicronsEndsThisFrame;
    
    save([p.tracksDir p.movieName '-skeletonData.mat'],...
        'allLengthsOfBacteriaInPixels','allLengthsOfBacteriaInMicrons',...
        'allPixelAreaOfBacterium','allSkeletonXYpoleToPole',...
        'allMinX','allMinY',...
        'allEdges','alldistanceAlongSkeletonPixels',...
        'allextrapolatedDistanceEndsPixels','allextrapolatedDistanceEndsMicrons');    
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
% figure(60); clf;
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