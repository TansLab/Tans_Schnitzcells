
function MW_straightenbacteria(p, frameRange, fluorColor) 

% MW simple straightening algorithm
% Thanks to Nicola Gritti for help with algorithm.
% 2016/04
%
% input:
%   p               our standard parameter struct
%   framerange      range of frames you want to process
%   fluorColor      one-letter abbreviation for used fluor: 'g','y','r','c'
%   extraOutput     optional to give more info (plots) to user
% Hard-coded parameters:
%   USEEXTENDEDSKELETON 
%                   whether to use skeleton+extrapolated ends, or only skeleton itself
%   SMOOTHELEMENTS  how many elements to take +/- each side for smoothing
% 
% output
%   - outputs to matlab file 
%       [outputDir p.movieDate p.movieName '_straightFluorData.mat']
%   - outputs plots to 
%       [outputDir p.movieDate p.movieName '_straightenedPlot_' fluorColor '_fr' sprintf('%03d',framenr) 'cell' sprintf('%03d',cellno)];
%
% parameters exported to .mat file
%   -  allskeletonDistance{framenr}{cellno} 
%         excluding extrapolated start
%   -  allmeanY{framenr}{cellno}          
%   -  allpeakIndices{framenr}{cellno}      
%   -  allpeakXPixels{framenr}{cellno}     
%         including extrapolated start
%   -  allpeakXMicrons{framenr}{cellno}    
%         including extrapolated start
%   -  allpeakmeanY{framenr}{cellno}        
%
% TODO
% Currently, bacteria is only straightened for branchless skeleton, not the
% extrapolated parts of the skeleton. Also peak distances are determined
% from this, but corrected using the extrapolated length of the extra ends.
%
% BACTERIAL WIDTH
% In principle it might be possible to determine the bacterial width rather
% precisely based on the segmentation files and the skeleton. I had
% introduced some code to address this, but I overlooked the facts that
% a) There is already a preliminary estimate of the bacterial width
% (area/skeleton lenght) which is used during straightening; as such the
% straightened skeleton might be biased by this estimate.
% b) If one estimates the width from the tangential lines in combination
% with the segmentation, then there should be a correction for the actual
% distance between pixels (i.e. the total number of pixels does not
% translate to a distance because they might be under an angle).
% Instead of addressing these issues, for now it is easier to go with the
% estimate (area/skeleton lenght) in our analyses. Bacterial width could be
% calculated more precisely later by extending this script.
% -MW

%% Parameters
% SMOOTHELEMENTS has been set to 0 since it is now already performed in
% NDL_lengthforfillamentedcells (was 8 before)
SMOOTHELEMENTS=8; % how many elements to take +/- each side 
USEEXTENDEDSKELETON = 1; % whether to use skeleton+extrapolated ends, or only skeleton itself

% Non-user parameters
if isfield(p,'extraOutput')
    extraOutput = p.extraOutput;
else
    extraOutput = 0;
end

outputDir = [p.analysisDir 'straightenedCells\'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if extraOutput   
    outputDirMore = [outputDir '\moreplots\'];
    if ~exist(outputDirMore,'dir')
        mkdir(outputDirMore);
    end
end

%% 
disp('Starting MW_straightenbacteria..');

%% Load skeleton data
load ([p.tracksDir p.movieName '-skeletonData.mat']);
if USEEXTENDEDSKELETON
    theSkeletons = allExtendedSkeletons;
else
    theSkeletons = allSkeletonXYpoleToPole; % allsmoothSkeletonXYpoleToPole
end

%% Prepare output parameters
allskeletonDistance         = {};
allrawFluor                 = {};
allbacterialWidthPixels           = {};
fluorInterpolatedValues     = {};
allpeakXPixels              = {};
allpeakXMicrons             = {};
allpeakX                    = {};
allpeakY                    = {};
allpeakIndices              = {};

%%
if extraOutput
    h104=figure();
else
    h104=figure('Visible','off');
    %set(h104,'Visible','off');
end

if extraOutput
    h105=figure();
else
    h105=figure('Visible','off');
    %set(h105,'Visible','off');
end

%% Loop over all frames
lastFrame=frameRange(end);
for framenr = frameRange
    
    disp(['Analyzing frame ' num2str(framenr) ' (highest framenr =' num2str(lastFrame) ').']);
    
    %% % Load data for current frame of the dataset
      
    %Load (corrected) fluorescence image
    fluorFileName = [p.tracksDir p.movieName 'Fluor_' fluorColor '_' sprintf('%03d',framenr) '.mat'];
    if exist(fluorFileName,'file')
        load(fluorFileName);
    else
        disp(['Fluor image not found, skipping this frame (' fluorFileName ').']);
        continue;
    end
    fluorImg = eval([fluorColor 'reg5']); % note is already shifted
        % xreg5 is corrected (and shifted) fluor image.
    
    % Load segmentation
    load ([p.segmentationDir p.movieName 'seg' sprintf('%03d',framenr) '.mat']);%,'Lc');    
        % Important contents is Lc, which respectively hold the
        % checked segmentation, and the fluorescence image, 
    segImg = Lc>0;
    
    %% Prepare loop over cellno's in this frame
    
    % get unique cellnos
    nonZeroIndices = (Lc(:)>0);
    allCellnos = transpose(unique(Lc(nonZeroIndices)));
    highestCellno = max(allCellnos);
    % Prepare output paramers    
    skeletonDistanceThisFrame       = cell(1,highestCellno);
    rawFluorThisFrame               = cell(1,highestCellno);
    bacterialWidthThisFrame         = cell(1,highestCellno);
    bacterialWidthThisFrameMicrons  = cell(1,highestCellno);
    fluorInterpolatedValues         = cell(1,highestCellno);
    peakXPixelsThisFrame            = cell(1,highestCellno);
    peakXMicronsThisFrame           = cell(1,highestCellno);
    peakYThisFrame                  = cell(1,highestCellno);
    peakIndicesThisFrame            = cell(1,highestCellno);
    
    disp(['    > ' num2str(highestCellno) ' cell(s) in this frame']);
    
    %% create figure
    if extraOutput
        h101=figure(101); clf;
        subplot(1,2,1);
        imshow(fluorImg,[]); hold on;
        subplot(1,2,2);
        imshow(phsub,[]); hold on;
    end
    %% loop over cellno's in this frame
    for cellno = allCellnos

        %% calculate/load some parameters for this cell
        currentSkeletonXYpoleToPole = theSkeletons{framenr}{cellno};
        currentminX = allMinX{framenr}(cellno);
        currentminY = allMinY{framenr}(cellno);
        currentEdges = allEdges{framenr}{cellno};
        currentSpacingInSkeleton = MW_distancesbetweenlinesegments(currentSkeletonXYpoleToPole(:,1)',currentSkeletonXYpoleToPole(:,2)');
        currentDistanceAlongSkeleton = cumsum(currentSpacingInSkeleton); % alldistanceAlongSkeletonPixels{framenr}{cellno};        
        
        averageBacterialWidth = allPixelAreaOfBacterium{framenr}(cellno)/allLengthsOfBacteriaInPixels{framenr}(cellno);
        halfAverageBacterialWidth = averageBacterialWidth/2;
        intAverageBacterialWidth = uint16(averageBacterialWidth);       

        %% plot skeleton
        if extraOutput
            figure(102); clf; hold on; axis equal;
            plot(currentSkeletonXYpoleToPole(:,1), currentSkeletonXYpoleToPole(:,2),'x')
            grid on
            grid minor
            set(gca,'XMinorTick','on','YMinorTick','on')
            hA = gca;
            hA.XRuler.MinorTick = [min(currentSkeletonXYpoleToPole(:,1)):1:max(currentSkeletonXYpoleToPole(:,1))];
            hA.YRuler.MinorTick = [min(currentSkeletonXYpoleToPole(:,2)):1:max(currentSkeletonXYpoleToPole(:,2))];
        end

        %% now average the elements
        windowArray = [-SMOOTHELEMENTS:SMOOTHELEMENTS];
        %smoothSkeleton = NaN(size(currentSkeletonXYpoleToPole,1)-2*SMOOTHELEMENTS+1,2)
        smoothSkeleton = NaN(size(currentSkeletonXYpoleToPole,1),2);
        nrIndicesInSkelet = size(currentSkeletonXYpoleToPole,1);
        % 1st few elements
        for i = 1:SMOOTHELEMENTS
            %%
            plusminus = i-1;
            smoothSkeleton(i,1) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
            smoothSkeleton(i,2) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
            %plot(smoothSkeleton(i,1), smoothSkeleton(i,2),'<k','MarkerFaceColor','r','MarkerSize',10); disp('REMOVE THIS LINE!!')
        end
        % main piece of skeleton
        for i = SMOOTHELEMENTS+1:nrIndicesInSkelet-SMOOTHELEMENTS

            smoothSkeleton(i,1) = mean(currentSkeletonXYpoleToPole(i+windowArray,1));
            smoothSkeleton(i,2) = mean(currentSkeletonXYpoleToPole(i+windowArray,2));

        end
        % last few elements
        for i = nrIndicesInSkelet-SMOOTHELEMENTS+1:nrIndicesInSkelet
            plusminus = nrIndicesInSkelet-i;
            smoothSkeleton(i,1) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
            smoothSkeleton(i,2) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
            %plot(smoothSkeleton(i,1), smoothSkeleton(i,2),'^k','MarkerFaceColor','b'); disp('REMOVE THIS LINE!!')
        end
        
        if extraOutput
            figure(102); hold on; 
            plot(smoothSkeleton(:,1), smoothSkeleton(:,2),'o');
            axis equal;
        end

         %% plot the skeleton on the original image
        if extraOutput
            figure(101); 
            subplot(1,2,1); hold on; 
            plot(smoothSkeleton(:,2)+currentminX, smoothSkeleton(:,1)+currentminY,'.')
            subplot(1,2,2); hold on; 
            plot(smoothSkeleton(:,2)+currentminX, smoothSkeleton(:,1)+currentminY,'.')
        end
        
        %% Now determine box coordinates and angles for all points

        pointsA = NaN(nrIndicesInSkelet-1,2);
        pointsB = NaN(nrIndicesInSkelet-1,2);
        tangentialLineXY1 = NaN(nrIndicesInSkelet-1,2);
        tangentialLineXY2 = NaN(nrIndicesInSkelet-1,2);
        tangentialLinesAllCoordsX = NaN(nrIndicesInSkelet-1,intAverageBacterialWidth);
        tangentialLinesAllCoordsY = NaN(nrIndicesInSkelet-1,intAverageBacterialWidth);
        vectonext = NaN(nrIndicesInSkelet-1,2);
        angles = NaN(nrIndicesInSkelet-1,1);
        for i = 1:(nrIndicesInSkelet-1)

            % two consecutive points
            point1 = smoothSkeleton(i,:);
            point2 = smoothSkeleton(i+1,:);

            % x2-x1 and y2-y1
            sideLength1 = point1(1)-point2(1);
            sideLength2 = point1(2)-point2(2);

            % angle between points
            theangle = atan(sideLength2/sideLength1);

            % box corner w. respect to x,y
            deltax = -sin(theangle)*halfAverageBacterialWidth;
            deltay = cos(theangle)*halfAverageBacterialWidth;

            %pointsA(i,:) = point1+[deltax, deltay];
            %pointsB(i,:) = point2-[deltax, deltay];

            %angles(i) = theangle;

            %vectonext(i,:)=[sideLength1,sideLength2];

            centerXY = (point1+point2)./2;
            
            tangentialLineXY1(i,:) = centerXY+[deltax, deltay];%+.5*[sideLength1 sideLength2];
            tangentialLineXY2(i,:) = centerXY-[deltax, deltay];%+.5*[sideLength1 sideLength2];

            tangentialLinesAllCoordsX(i,:) = uint16(round(linspace(tangentialLineXY1(i,1),tangentialLineXY2(i,1),intAverageBacterialWidth)));
            tangentialLinesAllCoordsY(i,:) = uint16(round(linspace(tangentialLineXY1(i,2),tangentialLineXY2(i,2),intAverageBacterialWidth)));

        end
        
        if extraOutput
            figure(102); hold on; 
            for i = 1:(nrIndicesInSkelet-1)
                %plot([smoothSkeleton(i,1),pointsA(i,1)], [smoothSkeleton(i,2), pointsA(i,2)],'-')
                %plot([smoothSkeleton(i+1,1),pointsB(i,1)], [smoothSkeleton(i+1,2), pointsB(i,2)],'-')

                plot([tangentialLineXY1(i,1),tangentialLineXY2(i,1)], [tangentialLineXY1(i,2),tangentialLineXY2(i,2)],'-')
                plot(tangentialLinesAllCoordsX(i,:),tangentialLinesAllCoordsY(i,:),'.')
                %rectangle('Position',[pointsA(i,1) pointsA(i,2) pointsB(i,1) pointsB(i,2)])
            end

            plot(currentEdges(:,1),currentEdges(:,2),'-');
        end

        %% Create a figure with the tangential lines on top of original
        if extraOutput
            figure(103); clf; hold on; 
            subplottight(2,1,1);
            imshow(phsub',[]); hold on;
            plot(currentSkeletonXYpoleToPole(:,1)+currentminY,currentSkeletonXYpoleToPole(:,2)+currentminX,'.')
            for i = 1:(nrIndicesInSkelet-1)
                plot(tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX,'.');
            end
            %
            subplottight(2,1,2);
            imshow(fluorImg',[]); hold on;
            plot(currentSkeletonXYpoleToPole(:,1)+currentminY,currentSkeletonXYpoleToPole(:,2)+currentminX,'.')
            for i = 1:(nrIndicesInSkelet-1)
                plot(tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX,'.');
            end
            
        end

        %% Now create the straightened bacteria
        % CAREFUL: the axis along the bacterium will be distorted! (Since the
        % distance between points on the skeleton is not equal.)

        straightenedBacterium       = NaN(intAverageBacterialWidth,nrIndicesInSkelet-1);
        straightenedBacteriumFluor  = NaN(intAverageBacterialWidth,nrIndicesInSkelet-1);
        straightenedBacteriumSeg    = NaN(intAverageBacterialWidth,nrIndicesInSkelet-1);
        theSizeBacteriaImage        = size(phsub); % size ph = size fluor        
        
        for i = 1:(nrIndicesInSkelet-1)
            %straightenedBacterium(:,i) = phsub(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+minX,tangentialLinesAllCoordsY(i,:)+minY));
            straightenedBacterium(:,i)      = phsub   (sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX));
            straightenedBacteriumFluor(:,i) = fluorImg(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX));
            straightenedBacteriumSeg(:,i)   = segImg(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX));
        end

        %{
        if extraOutput
            figure(h104); clf;
        else
            h=figure(h104);
            set(h,'Visible','off'); clf;
        end        
        %}
        
        distanceTangentials = currentDistanceAlongSkeleton-currentSpacingInSkeleton./2;        
        fluorIntensity = mean(straightenedBacteriumFluor);
        %bacterialWidth = sum(straightenedBacteriumSeg); 
            % parts referring to bacterial width have been removed, see
            % remark at top of this function
        
        %% Now create the distance-corrected straightened bacterium
        newXaxis = [ceil(distanceTangentials(1))+1:floor(distanceTangentials(end))-1];
        interpolatedFluorValues = NaN(1,numel(newXaxis));
        interpolatedWidthValues = NaN(1,numel(newXaxis));
        for i = 1:numel(newXaxis)
            x = newXaxis(i);
            % first element to the left of this x            
            leftIndexAll = find(distanceTangentials<x);
            leftIndex = leftIndexAll(end);
            % first element on the right of this x
            rightIndexAll = find(distanceTangentials>x); 
            rightIndex = rightIndexAll(1); % should be leftIndex+1
            % Determine the interpolated values
            DY = fluorIntensity(rightIndex)-fluorIntensity(leftIndex);
            DX = distanceTangentials(rightIndex)-distanceTangentials(leftIndex);
            dxprime = x-distanceTangentials(leftIndex);
            interpolatedFluorValues(i) = fluorIntensity(leftIndex)+(DY/DX)*dxprime;
            % Determine interpolated value for width
            %DY = bacterialWidth(rightIndex)-bacterialWidth(leftIndex);
            %interpolatedWidthValues(i) = bacterialWidth(leftIndex)+(DY/DX)*dxprime;
        end
        
        % Calculate peak parameters        
        [peakindexes, peakys] = peakfinder(interpolatedFluorValues);
        peakxs = newXaxis(peakindexes);
        
        %% Create plot
        % ===
        set(0, 'currentfigure', h104); clf;
        
        subplot(3,1,1);
        imshow(straightenedBacterium,[]);
        %===
        subplot(3,1,2);
        imshow(straightenedBacteriumFluor,[]);
        %===
        %{
        subplot(5,1,3);
        plot(currentDistanceAlongSkeleton','.')
        ylim([min(currentDistanceAlongSkeleton),max(currentDistanceAlongSkeleton)]);
        %pixelToPixelDistance = currentdistanceAlongSkeleton(1:end-1)-currentdistanceAlongSkeleton(2:end);
        %plot(pixelToPixelDistance','.')
        %ylim([0,2])
        xlim([0,numel(currentDistanceAlongSkeleton)]);
        ylabel('distance in pixels');
        MW_makeplotlookbetter(20);
        % === fluor intensity
        subplot(5,1,4);        
        plot(fluorIntensity','.','LineWidth',3); hold on;
        xlim([0,numel(currentDistanceAlongSkeleton)]);
        %ylim([min(currentdistanceAlongSkeleton),max(currentdistanceAlongSkeleton)]);
        ylim([0,max(fluorIntensity)*1.1])
        xlabel('pixels');
        ylabel('mean fluor (a.u.)');        
        plot(peakindexes, peakys,'o','MarkerSize',10,'LineWidth',3);
        MW_makeplotlookbetter(20);
        %}
        % == distance-corrected 
        subplot(3,1,3);
        plot(newXaxis,interpolatedFluorValues,'.','LineWidth',3); hold on;
        xlim([0,max(currentDistanceAlongSkeleton)]);
        %ylim([min(currentdistanceAlongSkeleton),max(currentdistanceAlongSkeleton)]);
        if any(interpolatedFluorValues>0)
            ylim([0,max(interpolatedFluorValues)*1.1])        
        end
        xlabel('pixels'); ylabel('mean fluor (a.u.)');      
        plot(peakxs, peakys,'o','MarkerSize',10,'LineWidth',3);
        %plot(peakindexes, peakys,'o','MarkerSize',10,'LineWidth',3);
        MW_makeplotlookbetter(20);
        % save plot        
        saveLocation = [outputDir p.movieDate p.movieName '_straightenedPlot_' fluorColor '_fr' sprintf('%03d',framenr) 'cell' sprintf('%03d',cellno)];
        saveas(gca, [saveLocation '.tif']);
    
        %% Figure of straightened segmentation for determining bacterial width
        %{
        set(0, 'currentfigure', h105); clf;
        
        subplot(4,1,1);
        imshow(straightenedBacterium,[]);
        subplot(4,1,2);
        imshow(straightenedBacteriumSeg, []);
        subplot(4,1,3)                 
        plot(bacterialWidth,'-','LineWidth',2);
        xlabel('Bacterial axis (px)'); ylabel('Width (px)');     
        subplot(4,1,4)
        plot(newXaxis.*p.micronsPerPixel,interpolatedWidthValues.*p.micronsPerPixel,'-','LineWidth',2);                        
        ylim([0,1]);
        xlabel('Bacterial axis (?m)'); ylabel('Width (?m)');      
        MW_makeplotlookbetter(12);
        saveLocation = [outputDir p.movieDate p.movieName '_straightenedSegPlot_fr' sprintf('%03d',framenr) 'cell' sprintf('%03d',cellno)];
        saveas(h105, [saveLocation '.tif']);
        %}
        %% /
    
        skeletonDistanceThisFrame{cellno}      = currentDistanceAlongSkeleton;
        rawFluorThisFrame{cellno}              = fluorIntensity;
        %bacterialWidthThisFrame{cellno}        = bacterialWidth;
        %bacterialWidthThisFrameMicrons{cellno} = bacterialWidth.*p.micronsPerPixel;        
        fluorInterpolatedValues{cellno}        = interpolatedFluorValues;
        peakIndicesThisFrame{cellno}           = peakindexes;
        peakXPixelsThisFrame{cellno}           = peakxs;
        peakXMicronsThisFrame{cellno}          = peakxs*p.micronsPerPixel;
        peakYThisFrame{cellno}                 = peakys;
        
    end
     
    allskeletonDistance{framenr} = skeletonDistanceThisFrame;
    allrawFluor{framenr}            = rawFluorThisFrame;
    %allbacterialWidthPixels{framenr}   = bacterialWidthThisFrame;
    %allbacterialWidthMicrons{framenr}   = bacterialWidthThisFrameMicrons;
    allFluorInterpolatedValues{framenr}   = fluorInterpolatedValues;
    allpeakIndices{framenr}      = peakIndicesThisFrame;
    allpeakXPixels{framenr}      = peakXPixelsThisFrame;
    allpeakXMicrons{framenr}     = peakXMicronsThisFrame;
    allpeakY{framenr}        = peakYThisFrame;    
    
    % save some figures if extra output was desired
    if extraOutput   
        saveas(h101, [outputDirMore 'skeletons_' num2str(framenr) '.tif']);
        saveas(h101, [outputDirMore 'skeletons_' num2str(framenr) '.fig']);
    end
    
end

%% Save data for later processing
saveLocationMatFile = [outputDir p.movieDate p.movieName '_straightFluorData.mat'];
    save(saveLocationMatFile,...
        'allskeletonDistance', 'allrawFluor', 'allFluorInterpolatedValues', 'allpeakXPixels',...
        'allpeakXMicrons', 'allpeakY');
    %, 'allbacterialWidthPixels', 'allbacterialWidthMicrons');   

disp(['Completely done straightening. Data saved to ' saveLocationMatFile]);    
    
end




%% 
%{
figure(); clf;
imshow(greg,[]);
%}

%% plot some boxes
%for i = 1:(nrIndicesInSkelet-1)
%{
for i = 50
    
    currentbox =    [ pointsA(i,:); ... %x1, y1
                      pointsA(i,:)+vectonext(i,:); ... %x2, y2
                      pointsB(i,:); ...%x3, y3
                      pointsB(i,:)-vectonext(i,:)]; %x4, y4
    
              for j = 1:3
                  line([currentbox(j,1),currentbox(j+1,1)],[currentbox(j,2),currentbox(j+1,2)],'Color',[1 1 1])
              end
              line([currentbox(1,1),currentbox(end,1)],[currentbox(1,2),currentbox(end,2)],'Color',[1 1 1])
end
%}

%% create rotated boxes
%{
% set box size again to test
warning('overwriting parameters that are set earlier');
BOXHEIGHT = 100;
BOXWIDTH = 10;


% create box
standardBox = ones(BOXHEIGHT,BOXWIDTH);
standardBox = padarray(standardBox,[0,(BOXHEIGHT-BOXWIDTH)/2],'both'); % make total square
% pad to allow for rotation
thePadding = round(sqrt(BOXHEIGHT^2+BOXWIDTH^2)-BOXHEIGHT);
standardBox = padarray(standardBox,[0,(BOXHEIGHT-BOXWIDTH)/2],'both');

for i=1:12
    theAngle = 360*i/12;
    rotatedBox = imrotate(standardBox,theAngle,'crop')
    figure(), imshow(rotatedBox,[]);
end
%}