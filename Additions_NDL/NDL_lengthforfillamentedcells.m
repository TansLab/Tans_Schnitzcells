

function [p,schnitzcells] = NDL_lengthforfillamentedcells(p) 


%% 
load('F:\Datasets\2016-03-11\pos4crop\segmentation\pos4cropseg191.mat')
load('F:\Datasets\2016-03-11\pos1crop\segmentation\pos1cropseg006.mat')
%load('F:\Datasets\2016-03-11\pos4crop\segmentation\pos4cropseg191.mat')

%% 
% Plots segmented image
figure(1); clf; imshow(Lc,[]);
cellnum=1;
[y,x] = find(Lc == cellnum);

%%
% Conversion
figure(2); clf;
plot(x,y,'.');
%%
% Makes image binary
zer=zeros(max(y),max(x));
for i=1:length(x)
    zer(y(i),x(i))=1;
end
bin_im=zer;
figure(3)
imshow(bin_im)
%%
% Skeletonizes image
BW = bwmorph(bin_im,'skel',Inf);
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