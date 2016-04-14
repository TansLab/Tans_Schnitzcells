

function [p,schnitzcells] = NDL_lengthforfillamentedcells(p) 

% function parameters set by user
micronsPerPixel=0.0431; % TODO! make this p.micronsPerPixel
AVERAGEBACTERIAWIDTH = .5;

% paramters calculated based on user-supplied parameters
averageBacterialWidthInPixel= AVERAGEBACTERIAWIDTH/micronsPerPixel;
paddingsize = round(averageBacterialWidthInPixel*4);

%% 
%load('F:\Datasets\2016-03-11\pos4crop\segmentation\pos4cropseg191.mat')
%load('F:\Datasets\2016-03-11\pos1crop\segmentation\pos1cropseg006.mat')
load('F:\Datasets\2016-03-11\pos4crop\segmentation\pos4cropseg191.mat')
% load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg233.mat'
% load 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-03-23\pos4crop\segmentation\pos4cropseg337.mat'


%% 
% Plots segmented image
cellnum=1;
figure(1); clf; imshow(Lc,[]);
[y,x] = find(Lc == cellnum);
% binary_image =  (Lc == cellnum); % immediate obtain binary, but no
% selection

%%
% Conversion
figure(2); clf;
plot(x,y,'.');

%%
% Select ROI and make image binary

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

% add padding to image (to avoid filters "seeing" edges)
bin_im=padarray(bin_im,[paddingsize,paddingsize]);

figure(3)
imshow(bin_im)

%% smooth using circular neighborhood
% Seems to make it a little better for my cells
SCALING = 1; % scaling is not so useful actually, 1 means no scaling is applied

work_img = bin_im;

% rescaling, maybe higher resolution better?
work_img = imresize(work_img,SCALING,'nearest');

% set circular neighborhood
% note that the size of the disk is pretty critical, and due to the size 
% of the bacterium is too small to have an actual effect
sizeOfDisk=round(averageBacterialWidthInPixel)*.5*SCALING;
se=strel('disk',sizeOfDisk);

% apply filters
%work_img = imerode(work_img,se);
%work_img = imdilate(work_img,se);
work_img = imclose(work_img,se); % works best

% show results
figure(); 
subplot(1,2,1);
imshow(bin_im);
subplot(1,2,2)
imshow(work_img);

if 0 % to turn on/off this filter
    bin_im = work_img;
end

%%
% Skeletonizes image
BW = bwmorph(bin_im,'skel',Inf);
%BW = voronoiSkel(bin_im); % downloaded this, but to tricky to get to work
%BW = skeleton(bin_im); % downloaded this, but to tricky to get to work

figure(4); clf;
imshow(BW)
imshow((bin_im+BW)/2,[])

%%
% Finds endings
ends = bwmorph(BW,'endpoints');
figure(50); clf;
imshow(ends)
%%
% Calculate number of ends
num_ends=sum(sum(ends,2))
%%
% Removes side-branches
count=0;
while num_ends>2
    BW = bwmorph(BW,'spur');
    count=count+1;
    ends = bwmorph(BW,'endpoints');
    num_ends=sum(sum(ends,2));
end
count
num_ends
BW1=BW;
figure(5); clf;
imshow(BW1)
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
end