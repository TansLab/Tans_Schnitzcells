function [theout] = progthreshsubtilis(im0,mask,sizefactor,invertit)
     % progthresh:
% inputs = im, nsteps

  %%%%%%%%%%%% directories & filenames:

% ima_dir = 'd:\hypermelow\matlab\bactrack\';
% pics = dir([ima_dir, '*.mat']);
% ima_N = length(pics);
% p = 50;
% fname = [ima_dir pics(p).name];
% load(fname); ima = double(phsub);
% im = ima - min2(ima);

%%%%%%%%%%%% progrthresh routine (input = im)

if nargin < 3,
sizefactor = 1;
end;
im = im0;
im = imresize(im,sizefactor,'bilinear');

minArea = 75;
maxArea = 2500;
minMinorAxis = 8;
maxMinorAxis = 21;
minMajorAxis = 17;
maxMajorAxis = 250;
minSolidity = 0.7; % was .8
minEccentricity = 0.6; % was .8
minRatio = 1.1; % the ratio of length to width
nsteps = 20;
maxsteps = 100;
im = medfilt2(im,[3 3]);

if nargin<2,
    mask=1;
end;


%minim = min(min(im)); maxim = median(im(:)); range = maxim - minim;
[muhat,sigmahat] = normfit(im(:));
medianim = median(im(:));
%minim = medianim - 4*sigmahat; maxim = medianim - 0.5*sigmahat; range = maxim - minim;
s = sort(im(:));

minim = s(20);
maxim = s(end-20);
range = maxim - minim;
%minim = 450; maxim = 456; range = maxim - minim;
%stepsize = max(range/nsteps,1);
stepsize = max(5,range/maxsteps);
kept = zeros([size(im)]);

% figure(1);
% imshow(im,[]);

Nims = length(minim:stepsize:maxim);
progthreshout = uint8(0);
progthreshout(size(im,1),size(im,2),Nims)=0;

jj = 0;
for T = minim:stepsize:maxim,
    jj = jj + 1;
[minim T maxim]

    % L = bwlabel(im < T);
    % just threshold for now- don't recolor before eroding
    %     L = bwlabel(im < T);
if invertit,
         L = im < T;
else
         L = im > T;
end;

    % erode blobs in case they're joined cells
    Le = imerode(L,strel('diamond',2));
    % recolor connected blobs
    Lec = bwlabel(Le);
    % now expand back to original size keeping unique cell identities
    L = carefuldilate(Lec,strel('diamond',1),3);
    L(mask==0)=0;

if max2(L) > 5,
         T = T + 2 * stepsize;
end;

r = regionprops(L, 'Area','MajorAxisLength','MinorAxisLength','Eccentricity','Solidity','EulerNumber');

% these filter parameters are just a first guess:
f = find(([r.Area]>minArea) & ...
         ([r.Area]<maxArea) & ...
         ([r.MinorAxisLength]>minMinorAxis) & ...
         ([r.MinorAxisLength]<maxMinorAxis) & ...
         ([r.MajorAxisLength]>minMajorAxis) & ...
         ([r.MajorAxisLength]<maxMajorAxis) & ...
         ([r.Solidity]>minSolidity) & ...
         ([r.Eccentricity]>minEccentricity) & ...
         ([r.EulerNumber]==1) & ...
         ([r.MajorAxisLength]./[r.MinorAxisLength] > minRatio));

Lkeep = L.*ismember(L,f);

% possible thing to do here: imclose or open to get more out of thresholds
%imopen(strel(Lkeep,'disk',2))
kept = Lkeep>0; % this keeps only the objects that pass the filter
progthreshout(:,:,jj) = uint8(kept);
end;

s = sum(progthreshout,3);
[es,ethresh] = edge(s,'log',0);
ethresh
esf = imfill(es,'holes');
theout = esf;
theout = imresize(theout,size(im0),'bilinear');

% the object we want to keep is the one closest to the middle:
L = bwlabel(theout);
r = regionprops(L,'Centroid');
if length(r)>0,
for i = 1:length(r),
         d2(i) = (r(i).Centroid(1)-size(theout,1)/2)^2 + (r(i).Centroid(2)-size(theout,2)/2)^2;
end;
f = find(d2==min(d2));
theout = L==f;
else
disp('nothing found');
end;

% theout = sum(progthreshout,3)>0;
% theout = imresize(theout,size(im0),'bilinear');
% oldfig = gcf; figure(201); imshow(theout);
% figure(oldfig);

% out1 = kept>0;
% out1 = imresize(out1, size(im0));
