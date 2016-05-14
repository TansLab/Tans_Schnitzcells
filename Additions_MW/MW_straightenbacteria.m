
function MW_straightenbacteria(p, frameRange, fluorColor) 

% MW simple straightening algorithm
% Thanks to Nicola Gritti for help with algorithm.
% 2016/04
%
% input:
%   p               our standard parameter struct
%   framerange      range of frames you want to process
%   fluorColor      one-letter abbreviation for used fluor: 'g','y','r','c'
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

%% Parameters
SMOOTHELEMENTS=8; % how many elements to take +/- each side

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

%% 
disp('Starting MW_straightenbacteria..');

%% Load skeleton data
load ([p.tracksDir p.movieName '-skeletonData.mat']);

%% Prepare output parameters
allskeletonDistance = {};
allmeanY            = {};
allpeakXPixels      = {};
allpeakXMicrons     = {};
allpeakmeanY        = {};
allpeakIndices      = {};

%% Loop over all frames
lastFrame=frameRange(end);
for framenr = frameRange
    
    disp(['Analyzing frame ' num2str(framenr) ' (highest framenr =' num2str(lastFrame) ').']);
    
    %% % Load data for current frame of the dataset
    
    %e.g. load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg337.mat'
    % Load segmentation
    load ([p.segmentationDir p.movieName 'seg' sprintf('%03d',framenr) '.mat']);%,'Lc');
        % Important contents is Lc, which respectively hold the
        % checked segmentation, and the fluorescence image, 
    %Load (corrected) fluorescence image
    load([p.tracksDir p.movieName 'Fluor_' fluorColor '_' sprintf('%03d',framenr) '.mat'],[fluorColor 'reg5']);
    %G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-04-07\pos1crop\data
    fluorImg = eval([fluorColor 'reg5']); % note is already shifted
        % xreg5 is corrected (and shifted) fluor image.
        
    
    %% Loop over all cell numbers in this frame
    
    % get unique cellnos
    nonZeroIndices = (Lc(:)>0);
    allCellnos = transpose(unique(Lc(nonZeroIndices)));
    highestCellno = max(allCellnos);
    % Prepare output paramers    
    skeletonDistanceThisFrame = cell(1,highestCellno);
    meanYThisFrame            = cell(1,highestCellno);
    peakXPixelsThisFrame      = cell(1,highestCellno);
    peakXMicronsThisFrame     = cell(1,highestCellno);
    peakmeanYThisFrame        = cell(1,highestCellno);
    peakIndicesThisFrame      = cell(1,highestCellno);
    
    disp(['    > ' num2str(highestCellno) ' cell(s) in this frame']);
    
    % loop
    for cellno = allCellnos

        %% calculate/load some parameters for this cell
        currentSkeletonXYpoleToPole = allSkeletonXYpoleToPole{framenr}{cellno};
        currentminX = allMinX{framenr}(cellno);
        currentminY = allMinY{framenr}(cellno);
        currentarray2 = allarray2{framenr}{cellno};
        currentdistanceAlongSkeleton = alldistanceAlongSkeletonPixels{framenr}{cellno};
        
        averageBacterialWidth = allPixelAreaOfBacterium{framenr}(cellno)/allLengthsOfBacteriaInPixels{framenr}(cellno);
        halfAverageBacterialWidth = averageBacterialWidth/2;
        intAverageBacterialWidth = uint16(averageBacterialWidth);

        %% plot the skeleton on the original image
        if extraOutput
            figure(101); clf; hold on; 
            subplot(1,2,1);
            imshow(fluorImg,[]); hold on;
            plot(currentSkeletonXYpoleToPole(:,2)+currentminX, currentSkeletonXYpoleToPole(:,1)+currentminY,'.')
            subplot(1,2,2);
            imshow(phsub,[]); hold on;
            plot(currentSkeletonXYpoleToPole(:,2)+currentminX, currentSkeletonXYpoleToPole(:,1)+currentminY,'.')
        end

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
            plusminus = i-1;
            smoothSkeleton(i,1) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
            smoothSkeleton(i,2) = mean(currentSkeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
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
        end
        
        if extraOutput
            figure(102); hold on; 
            plot(smoothSkeleton(:,1), smoothSkeleton(:,2),'o')
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
            vec1 = smoothSkeleton(i,:);
            vec2 = smoothSkeleton(i+1,:);

            % x2-x1 and y2-y1
            sideLength1 = vec1(1)-vec2(1);
            sideLength2 = vec1(2)-vec2(2);

            % angle between points
            theangle = atan(sideLength2/sideLength1);

            % box corner w. respect to x,y
            deltax = -sin(theangle)*halfAverageBacterialWidth;
            deltay = cos(theangle)*halfAverageBacterialWidth;

            pointsA(i,:) = smoothSkeleton(i,:)+[deltax, deltay];
            pointsB(i,:) = smoothSkeleton(i+1,:)-[deltax, deltay];

            angles(i) = theangle;

            vectonext(i,:)=[sideLength1,sideLength2];

            tangentialLineXY1(i,:) = smoothSkeleton(i,:)+[deltax, deltay]+.5*[sideLength1 sideLength2];
            tangentialLineXY2(i,:) = smoothSkeleton(i,:)-[deltax, deltay]+.5*[sideLength1 sideLength2];

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

            plot(currentarray2(:,1),currentarray2(:,2),'-');
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
        theSizeBacteriaImage        = size(phsub); % size ph = size fluor        
        
        for i = 1:(nrIndicesInSkelet-1)
            %straightenedBacterium(:,i) = phsub(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+minX,tangentialLinesAllCoordsY(i,:)+minY));
            straightenedBacterium(:,i)      = phsub   (sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX));
            straightenedBacteriumFluor(:,i) = fluorImg(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+currentminY,tangentialLinesAllCoordsY(i,:)+currentminX));
        end

        if extraOutput
            figure(104); clf;
        else
            h=figure(104);
            set(h,'Visible','off'); clf;
        end

        % Calculate peak parameters
        fluorIntensity = mean(straightenedBacteriumFluor);
        [peakindexes, peakys] = peakfinder(fluorIntensity);
        
        % ===
        subplot(4,1,1);
        imshow(straightenedBacterium,[]);
        %===
        subplot(4,1,2);
        imshow(straightenedBacteriumFluor,[]);
        %===
        subplot(4,1,3);
        plot(currentdistanceAlongSkeleton','.')
        ylim([min(currentdistanceAlongSkeleton),max(currentdistanceAlongSkeleton)]);
        %pixelToPixelDistance = currentdistanceAlongSkeleton(1:end-1)-currentdistanceAlongSkeleton(2:end);
        %plot(pixelToPixelDistance','.')
        %ylim([0,2])
        xlim([0,numel(currentdistanceAlongSkeleton)]);
        ylabel('distance in pixels');
        MW_makeplotlookbetter(20);
        % === fluor intensity
        subplot(4,1,4);        
        plot(fluorIntensity','.','LineWidth',3); hold on;
        xlim([0,numel(currentdistanceAlongSkeleton)]);
        %ylim([min(currentdistanceAlongSkeleton),max(currentdistanceAlongSkeleton)]);
        ylim([0,max(fluorIntensity)*1.1])
        xlabel('pixels');
        ylabel('mean fluor (a.u.)');        
        plot(peakindexes, peakys,'o','MarkerSize',10,'LineWidth',3);
        MW_makeplotlookbetter(20);            
        saveLocation = [outputDir p.movieDate p.movieName '_straightenedPlot_' fluorColor '_fr' sprintf('%03d',framenr) 'cell' sprintf('%03d',cellno)];
        saveas(gca, [saveLocation '.tif']);
    
        %% /
    
        skeletonDistanceThisFrame{cellno} = currentdistanceAlongSkeleton;
        meanYThisFrame{cellno}            = fluorIntensity;
        peakIndicesThisFrame{cellno}      = peakindexes;
        peakXPixelsThisFrame{cellno}      = allextrapolatedDistanceEndsPixels{framenr}(1,cellno) + ...
                                                currentdistanceAlongSkeleton(peakindexes);
        peakXMicronsThisFrame{cellno}     = allextrapolatedDistanceEndsMicrons{framenr}(1,cellno) + ...
                                                currentdistanceAlongSkeleton(peakindexes)*p.micronsPerPixel;
        peakmeanYThisFrame{cellno}        = peakys;
        
    end
     
    allskeletonDistance{framenr} = skeletonDistanceThisFrame;
    allmeanY{framenr}            = meanYThisFrame;
    allpeakIndices{framenr}      = peakIndicesThisFrame;
    allpeakXPixels{framenr}      = peakXPixelsThisFrame;
    allpeakXMicrons{framenr}     = peakXMicronsThisFrame;
    allpeakmeanY{framenr}        = peakmeanYThisFrame;    
    
    
    
end

% Save data for later processing
saveLocationMatFile = [outputDir p.movieDate p.movieName '_straightFluorData.mat'];
    save(saveLocationMatFile,...
        'allskeletonDistance', 'allmeanY', 'allpeakXPixels',...
        'allpeakXMicrons', 'allpeakmeanY');

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