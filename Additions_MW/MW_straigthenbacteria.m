

% MW simple straightening algorithm

SMOOTHELEMENTS=8; % how many elements to take +/- each side
BOXHEIGHT = 10; % BOXHEIGHT>=BOXWIDTH and should be even nr.
BOXWIDTH = 2; % note that BOXHEIGHT-BOXWIDTH should be an even nr

averageBacterialWidth = pixelAreaOfBacterium/lengthOfBacteriaInPixels;
halfAverageBacterialWidth = averageBacterialWidth/2;
intAverageBacterialWidth = uint16(averageBacterialWidth);

%% plot the skeleton on the original image
figure(101); clf; hold on; 
imshow(phsub,[]); hold on;
plot(skeletonXYpoleToPole(:,2)+minX, skeletonXYpoleToPole(:,1)+minY,'.')

%% plot skeleton
figure(102); clf; hold on; axis equal;
plot(skeletonXYpoleToPole(:,1), skeletonXYpoleToPole(:,2),'.')
grid on
grid minor
set(gca,'XMinorTick','on','YMinorTick','on')
hA = gca;
hA.XRuler.MinorTick = [min(skeletonXYpoleToPole(:,1)):1:max(skeletonXYpoleToPole(:,1))];
hA.YRuler.MinorTick = [min(skeletonXYpoleToPole(:,2)):1:max(skeletonXYpoleToPole(:,2))];

%% now average the elements
windowArray = [-SMOOTHELEMENTS:SMOOTHELEMENTS];
%smoothSkeleton = NaN(size(skeletonXYpoleToPole,1)-2*SMOOTHELEMENTS+1,2)
smoothSkeleton = NaN(size(skeletonXYpoleToPole,1),2);
nrIndicesInSkelet = size(skeletonXYpoleToPole,1);
% 1st few elements
for i = 1:SMOOTHELEMENTS
    plusminus = i-1;
    smoothSkeleton(i,1) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
    smoothSkeleton(i,2) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
end
% main piece of skeleton
for i = SMOOTHELEMENTS+1:nrIndicesInSkelet-SMOOTHELEMENTS
       
    smoothSkeleton(i,1) = mean(skeletonXYpoleToPole(i+windowArray,1));
    smoothSkeleton(i,2) = mean(skeletonXYpoleToPole(i+windowArray,2));
    
end
% last few elements
for i = nrIndicesInSkelet-SMOOTHELEMENTS+1:nrIndicesInSkelet
    plusminus = nrIndicesInSkelet-i;
    smoothSkeleton(i,1) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),1));
    smoothSkeleton(i,2) = mean(skeletonXYpoleToPole(i+windowArray(SMOOTHELEMENTS+1-plusminus:SMOOTHELEMENTS+1+plusminus),2));
end
figure(102); hold on; 
plot(smoothSkeleton(:,1), smoothSkeleton(:,2),'o')

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

figure(102); hold on; 
for i = 1:(nrIndicesInSkelet-1)
    %plot([smoothSkeleton(i,1),pointsA(i,1)], [smoothSkeleton(i,2), pointsA(i,2)],'-')
    %plot([smoothSkeleton(i+1,1),pointsB(i,1)], [smoothSkeleton(i+1,2), pointsB(i,2)],'-')

    plot([tangentialLineXY1(i,1),tangentialLineXY2(i,1)], [tangentialLineXY1(i,2),tangentialLineXY2(i,2)],'-')
    plot(tangentialLinesAllCoordsX(i,:),tangentialLinesAllCoordsY(i,:),'.')
    %rectangle('Position',[pointsA(i,1) pointsA(i,2) pointsB(i,1) pointsB(i,2)])
end

plot(array2(:,1),array2(:,2),'-');

%% Create a figure with the tangential lines on top of original
figure(103); clf; hold on; 
imshow(phsub',[]); hold on;
plot(skeletonXYpoleToPole(:,1)+minY,skeletonXYpoleToPole(:,2)+minX,'.')
for i = 1:(nrIndicesInSkelet-1)
    plot(tangentialLinesAllCoordsX(i,:)+minY,tangentialLinesAllCoordsY(i,:)+minX,'.');
end

%% Now create the straightened bacteria
% CAREFUL: the axis along the bacterium will be distorted! (Since the
% distance between points on the skeleton is not equal.)

straightenedBacterium = NaN(intAverageBacterialWidth,nrIndicesInSkelet-1);
theSizeBacteriaImage=size(phsub);

for i = 1:(nrIndicesInSkelet-1)
    %straightenedBacterium(:,i) = phsub(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+minX,tangentialLinesAllCoordsY(i,:)+minY));
    straightenedBacterium(:,i) = phsub(sub2ind(theSizeBacteriaImage,tangentialLinesAllCoordsX(i,:)+minY,tangentialLinesAllCoordsY(i,:)+minX));
end

figure(); clf;
imshow(straightenedBacterium,[]);




%% 
figure(); clf;
imshow(greg,[]);

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