function shadingImage = MW_shading_averaging(myFolder, myImages, NeighborAveraging)
% function shadingImage = MW_shading_averaging(myFolder, myImages, NeighborAveraging)
% MW shading images script
%
% Creates median image of supplied list of images, if NeighborAveraging =
% [X, X] with X>1 also smooths image.
% Creates also two subsets (even and uneven) and supplies graphs to compare
% those two images.
%
% Example list of input parameters and calling function:
% myFolder = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-05-07_shading_images\engfp\';
% myImages = { 'gfp20ms-r1.tif', 'gfp20ms-r2.tif',... };
% NeighborAveraging = 5;
% shadingImage = MW_shading_averaging(myFolder, myImages, NeighborAveraging)
%
% Output is img that can be used as shading image, shuold be saves as
% .mat file.
%%%%%%%%%%%%%%%%%%

% Parameter settings
HISTN=200; % number of bins used for histogram.

%% Loading and processing images

% Administration
nrImages = numel(myImages);
set1=[1:2:nrImages]; % Set of uneven images
set2=[2:2:nrImages]; % Set of even images

% Obtain dimention of image
firstImage = imread([myFolder myImages{1}]);
myDimension = size(firstImage);

% Create storage for all images
allImages = zeros(nrImages,myDimension(1), myDimension(2));
normAllImages = zeros(nrImages,myDimension(1), myDimension(2));

% Loop over all given images
for i = 1:nrImages
    % load image
    currentImage = double(imread([myFolder myImages{i}]));
    % store raw image
    allImages(i,:,:) = currentImage;
    % store filtered and normalized image    
    currentImage = currentImage/median(currentImage(:)); % normalizing by median works better than max
    currentImage = medfilt2(currentImage,[NeighborAveraging,NeighborAveraging]); % median filters decreases noise
    normAllImages(i,:,:) = currentImage;
end

%% plot all histograms
figure; clf; hold on;
for i = 1:nrImages
    [myHist, xHist] = hist(normAllImages(i,:),HISTN);
    plot(xHist, myHist, '-');
end


%% calculate medians

% median images (also for two subsets)
medianImageAll =  squeeze(median(normAllImages(:,:,:))); % squeeze removes singleton dim
medianImage1 = squeeze(median(normAllImages(set1,:,:))); % squeeze removes singleton dim
medianImage2 = squeeze(median(normAllImages(set2,:,:))); % squeeze removes singleton dim

% Create RGB version of those
RBGmedianImageAll = repmat(medianImageAll,[1 1 3]);  
RBGmedianImage1 = repmat(medianImage1,[1 1 3]);  
RBGmedianImage2 = repmat(medianImage2,[1 1 3]);  

% show histogram of the three created imgs
figure; clf; hold on;
[myMedianHist1, xMedianHist1] = hist(medianImage1(:),HISTN);
[myMedianHist2, xMedianHist2] = hist(medianImage2(:),HISTN);
[myMedianHistAll, xMedianHistAll] = hist(medianImageAll(:),HISTN);
plot(xMedianHist1,myMedianHist1,'b-');
plot(xMedianHist2,myMedianHist2,'r-');
plot(xMedianHist2,myMedianHistAll,'k-');
xlabel('Intensity');
ylabel('Count');
legend({'Uneven subset','Even subset','All imgs'},'Location','NorthWest')

%% show images
figure(); clf; hold on; 
subplottight(1,3,1);
imshow(RBGmedianImage1); title('based subset uneven');
subplottight(1,3,2);
imshow(RBGmedianImage2); title('even subset');
subplottight(1,3,3);
imshow(RBGmedianImageAll); title('median all');

% calculate percentage-wise-difference between subsets
deviationImage = medianImage1./medianImage2;

% create histogram of that and show user
figure; 
[myDevHist1, xDevHist1] = hist(deviationImage(:),HISTN);
plot(xDevHist1, myDevHist1);
xlim([0.8,1.2]);
xlabel('distribution (matrix subset1)/(matrix subset2) values');
ylabel('count');

%% Final output   
shadingImage = medianImageAll;

end
    
    




