

function [p,schnitzcells] = NDL_lengthforfillamentedcells(p) 

% function parameters set by user
micronsPerPixel=0.0431; % TODO! make this p.micronsPerPixel
AVERAGEBACTERIAWIDTH = .5;
EXTRAPOLATIONLENGTH = 30;

% parameters calculated based on user-supplied parameters
averageBacterialWidthInPixel= AVERAGEBACTERIAWIDTH/micronsPerPixel;
paddingsize = round(averageBacterialWidthInPixel*4);
extrapolationLength=round(averageBacterialWidthInPixel);

% if isfield(p,'extraoutput')
%     % Plot with outline of all bacteria and extended skeletons
%     plot() 
%     saveas([p.analysisDir '/lengthNick/' num2str(frame) num2str(cellno) '.tif'])
% end

%% % Loads one frame of a dataset
%load('F:\Datasets\2016-03-11\pos1crop\segmentation\pos1cropseg101.mat')
%load('F:\Datasets\2016-03-11\pos4crop\segmentation\pos4cropseg191.mat')
% load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg233.mat'
 load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg337.mat'

%% % Plots segmented image of particular cell
cellnum=1;
figure(1); clf; imshow(Lc,[]);
[y,x] = find(Lc == cellnum);
% binary_image =  (Lc == cellnum); % immediate obtain binary, but no
% selection

%% % Conversion
figure(2); clf;
plot(x,y,'.');

%% % Select ROI and make image binary
% administration required to select ROI 
minY = min(y);
minX = min(x);
sizeY = max(y)-minY+1;
sizeX = max(x)-minX+1;

% zeros
zer=zeros(sizeY,sizeX);

for i=1:length(x)
    zer(y(i)-minY+1,x(i)-minX+1)=1;
end

bin_im=zer;

% % add padding to image (to avoid filters "seeing" edges)
% bin_im=padarray(bin_im,[paddingsize,paddingsize]);

figure(3)
imshow(bin_im)
pixelAreaOfBacterium = sum(bin_im(:));

% %% smooth using circular neighborhood
% % Seems to make it a little better for my cells
% SCALING = 1; % scaling is not so useful actually, 1 means no scaling is applied
% 
% work_img = bin_im;
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
% imshow(bin_im);
% subplot(1,2,2)
% imshow(work_img);
% 
% if 0 % to turn on/off this filter
%     bin_im = work_img;
% end

%% % Skeletonizes image
BW = bwmorph(bin_im,'skel',Inf);
%BW = voronoiSkel(bin_im); % downloaded this, but to tricky to get to work
%BW = skeleton(bin_im); % downloaded this, but to tricky to get to work

figure(4); clf;
imshow(BW)
imshow((bin_im+BW)/2,[])
%% XXX
test=bwconncomp(BW)
%% % Finds edges
edges = bin_im-bwmorph(bin_im,'erode');
figure(); clf;
imshow(edges+BW)
%% % Finds endings
ends_before = bwmorph(BW,'endpoints');
figure(50); clf;
imshow(ends_before)
%% % Calculate number of ends
num_ends_before=sum(sum(ends_before,2))
%% Finds x & y values of endings before removing branches
k=1;
xyends_before=zeros(1,2,num_ends_before);
for i=1:sizeX
    for j=1:sizeY
        if ends_before(j,i)==1
            xyends_before(:,:,k)=[j,i];
            k=k+1;
        end
    end
end
xyends_before
%% % Disconnects at branch points
disc = BW-bwmorph(BW,'branchpoints');
figure(); clf;
imshow(disc)
%% % XXX
[xx,yy]=find(BW==1);
plot(xx,yy,'.')
func=csaps(xx,yy);
extra=fnxtr(func)
figure()
fnplt(extra)
% BWfit=fit(xx,yy,'poly9')
% plot(BWfit)
% BWspline=spline(xx,yy)
%% % Removes side-branches
count=0;
num_ends=num_ends_before;
while num_ends>2
    BW = bwmorph(BW,'spur');
    count=count+1;
    ends = bwmorph(BW,'endpoints');
    num_ends=sum(sum(ends,2));
end
BW1=BW;
count
num_ends
%% % Finds endings
ends = bwmorph(BW1,'endpoints');
figure(51); clf;
imshow(ends)
%% Finds x & y values of endings after removing branches
l=1;
xyends=zeros(1,2,num_ends);
for i=1:sizeX
    for j=1:sizeY
        if ends(j,i)==1
            xyends(:,:,l)=[j,i];
            l=l+1;
        end
    end
end
xyends
%% Gets x & y values of the branchless skeleton, and sets them in an array

% Get "boundaries" of the line (result should make a loop around the line,
% but from an arbitrary point)
skeletonBoundary=bwboundaries(BW1,8); 
skeletonBoundary=skeletonBoundary{1,1};

% paste this loop 2x behind itself
twotimesskeletonBoundary = [skeletonBoundary; skeletonBoundary];
leftend = xyends(:,:,1);

% Now find one of the bacterial poles
poleIndex = find(skeletonBoundary(:,1)==leftend(1) & skeletonBoundary(:,2)==leftend(2));

% Now get skeletonXYpoleToPole(i,:)=v(i), with v(i,:)=[x(i),y(i)],
% point i along the skeleton
skeletonXYpoleToPole=twotimesskeletonBoundary(poleIndex:poleIndex+round((size(skeletonBoundary,1)+1)/2),:);

% for i=1:sizeX
%     for i=1:sizeY
%         
% if ends==1
%     xends=
figure(5); clf;
imshow((BW1+bin_im)/2)
%% % Gets x & y values of the segmented edges, and plots them
edges2=bwboundaries(edges,8);
array2=edges2{1,1};
figure(71)
plot(array2(:,1),array2(:,2))
hold on
plot(skeletonXYpoleToPole(:,1),skeletonXYpoleToPole(:,2))
%% % Sets length of dataset (coming from branchless skeleton) which gets extrapolated
% extra_length=30;
% if length(array)<extra_length
%     extra_length=length(array)
% end
extra_length = min(EXTRAPOLATIONLENGTH, length(skeletonXYpoleToPole));
% vq3 = interp1(array(1:20,1),array(1:20,2),'pchip');
% bla=bspline(array(1:50,1),array(1:50,2));
%% % Extrapolates one end
func=csaps(skeletonXYpoleToPole(1:extra_length,1),skeletonXYpoleToPole(1:extra_length,2)); % TODO MAYBE USE OTHER (POLY)FIT?
extra=fnxtr(func,2);
figure()
extrapolatedSkeleton1 = fnplt(extra,[skeletonXYpoleToPole(1,1)-(count+extrapolationLength) skeletonXYpoleToPole(1,1)+(count+extrapolationLength)]).'
fnplt(extra,[skeletonXYpoleToPole(1,1)-(count+extrapolationLength) skeletonXYpoleToPole(1,1)+(count+extrapolationLength)])
%% % Extrapolates other end
func2=csaps(skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extra_length-1):length(skeletonXYpoleToPole),1),skeletonXYpoleToPole(length(skeletonXYpoleToPole)-(extra_length-1):length(skeletonXYpoleToPole),2));
extrapolatedSkeleton2=fnxtr(func2);
figure()
points2 = fnplt(extrapolatedSkeleton2,[skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)-(count+extrapolationLength) skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)+(count+extrapolationLength)]).'
fnplt(extrapolatedSkeleton2,[skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)-(count+extrapolationLength) skeletonXYpoleToPole(length(skeletonXYpoleToPole),1)+(count+extrapolationLength)])
%% % Plot extrapolations and segmented edges
figure(72)
plot(array2(:,1),array2(:,2))
hold on
plot(extrapolatedSkeleton1(:,1),extrapolatedSkeleton1(:,2))
hold on
plot(points2(:,1),points2(:,2))
hold on
plot(skeletonXYpoleToPole(:,1),skeletonXYpoleToPole(:,2))
%% % Determine intersection point and with that the correction length for one end
disx=zeros(length(extrapolatedSkeleton1),length(array2));
disy=zeros(length(extrapolatedSkeleton1),length(array2));
distot=zeros(length(extrapolatedSkeleton1),length(array2));
for i=1:length(extrapolatedSkeleton1)
    for j=1:length(array2)
        disx(i,j)=abs(array2(j,1)-extrapolatedSkeleton1(i,1));
        disy(i,j)=abs(array2(j,2)-extrapolatedSkeleton1(i,2));
        distot(i,j)=disx(i,j)+disy(i,j);
    end
end
mini=min(min(distot))
for i=1:length(extrapolatedSkeleton1)
    for j=1:length(array2)
        if distot(i,j)==mini
            icor=i;
            jcor=j;
        end
    end
end
array2(jcor,:)
extrapolated_intersection=extrapolatedSkeleton1(icor,:)
xyends(1,:,1)
extra_dist11=pdist2(extrapolated_intersection,xyends(1,:,1));
extra_dist12=pdist2(extrapolated_intersection,xyends(1,:,num_ends));
extra_dist1=min([extra_dist11 extra_dist12])
%% % Determine other intersection point and with that the correction length for the other end
disx2=zeros(length(points2),length(array2));
disy2=zeros(length(points2),length(array2));
distot2=zeros(length(points2),length(array2));
for i=1:length(points2)
    for j=1:length(array2)
        disx2(i,j)=abs(array2(j,1)-points2(i,1));
        disy2(i,j)=abs(array2(j,2)-points2(i,2));
        distot2(i,j)=disx2(i,j)+disy2(i,j);
    end
end
mini2=min(min(distot2))
for i=1:length(points2)
    for j=1:length(array2)
        if distot2(i,j)==mini2
            icor2=i;
            jcor2=j;
        end
    end
end
array2(jcor2,:)
extrapolated_intersection2=points2(icor2,:)
xyends(1,:,num_ends)
extra_dist21=pdist2(extrapolated_intersection2,xyends(1,:,1));
extra_dist22=pdist2(extrapolated_intersection2,xyends(1,:,num_ends));
extra_dist2=min([extra_dist21 extra_dist22])
%% Calculates length of branchless skeleton, and the total estimated length (by adding the extrapolated lengths of the ends)
distance_mask=ends;
extract_end=xyends(1,:,1);
distance_mask(extract_end(1),extract_end(2))=0;
D=bwdistgeodesic(BW1,distance_mask,'quasi-euclidean');
dist_BW1=max(max(D))
lengthOfBacteriaInPixels=dist_BW1+extra_dist1+extra_dist2;
% int2=zeros(1,length(points));
% for i=1:length(points)
%     int2(1,i)=min(disx(i,:))+min(disy(i,:));
% end
% int2;
% %
% Alternative script to get rid of side-branches
% skel= bwmorph(bin_im,'skel',Inf);
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

end