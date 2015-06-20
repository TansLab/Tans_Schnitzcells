function outputImage = MW_background_averaging(myFolder, myImages, settings)
% function outputImage = MW_background_averaging(myFolder, myImages, settings)
% MW shading images script
%
% Creates median image of supplied list of images, if NeighborAveraging =
% [X, X] with X>1 also smooths image.
% Creates also two subsets (even and uneven) and supplies graphs to compare
% those two images.
%
% Input parameters 
% - myFolder            Folder where to find images
% - myImages            List of image filenames. Set to 'auto' to search for
%                       images.
% - settings            Struct for optional settings. Default values
%                       hard-coded.
%   settings.histN                  Number of bins used for histogram
%   settings.maxImages              Matlab doesn't handle loading a large nr of images
%                                   very well. So this is capped!
%   settings.neighborAveraging      Sets degree of smoothing. (How many neighbors to
%                                   take into account for neighborhood averaging of 
%                                   pixels.) Only use for shading, not for
%                                   flatfield!
%   settings.useNormalizedImages    Whether imgs should be normalized.
%   settings.backgroundImage        If a (same sized) background image
%                                   (sometimes called 'flatfield' throughout the script)
%                                   is given, this will be substracted from
%                                   the input images.
%
% Example list of input parameters and calling function:
% %%% Background image
%   myFolder = 'F:\A_Tans0_step1_incoming_not_backed_up\2015-05-27_background\OPT-flatfield\10ms\'
%   myImages = 'auto'; 
%   settings.neighborAveraging = 0;
%   settings.useNormalizedImages = 0;
%   settings.maxImages = 500;
%   settings.histN = 50;
%   outputImage = MW_background_averaging(myFolder, myImages, settings);
%
%
% Output is img that can be used as shading image, shuold be saves as
% .mat file.
%%%%%%%%%%%%%%%%%%

% Parameter settings
if ~isfield(settings, 'useNormalizedImages')
    settings.useNormalizedImages = 1; % whether normalized imgs should be used
        % MW TODO OPTIMALIZATION (optional) note that code creates 
        % normalized version nevertheless
    warning('useNormalizedImages set to default value 1');
end
if ~isfield(settings, 'histN')
    settings.histN = 50; % number of bins used for histogram.
    warning('histN set to default value 50');
end
if ~isfield(settings, 'maxImages')
    settings.maxImages = 200; % max. number of images that are used
    warning('maxImages set to default value 200');    
end
if ~isfield(settings, 'neighborAveraging')
    settings.neighborAveraging = 5; % max. number of images that are used
    warning('neighborAveraging set to default value 5');    
end 


%% Loading and processing images

% Automatically create list images if desired
if ~iscell(myImages), if (myImages == 'auto')
    myImages = dir(myFolder);     
    myImages = {myImages(3:end).name}; % get names and filter out '.' and '..'
end, end

% cap number of images used
if numel(myImages) > settings.maxImages
    myImages = myImages(1:settings.maxImages); 
end; 

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
    % substract background if given
    if isfield(settings, 'backgroundImage')
        currentImage = currentImage-settings.backgroundImage;
        disp('Substracting background.');
    end
    % store raw image
    allImages(i,:,:) = currentImage;
    % store filtered and normalized image    
    currentImage = currentImage/median(currentImage(:)); % normalizing by median works better than max    
    if settings.neighborAveraging>0 % 
        % smooth
        currentImage = medfilt2(currentImage,[settings.neighborAveraging,settings.neighborAveraging]); 
    end
    normAllImages(i,:,:) = currentImage;
    % update user
    disp(['Loading ' num2str(i) '/' num2str(nrImages) '.']);
end

%% plot all histograms
figure; clf; 
subplot(2,1,1); hold on;
xlabel('Intensity distribution raw imgs');
ylabel('Count');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);
subplot(2,1,2); hold on;
xlabel('Normalized, processed');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);
% Create colormap
myColors = distinguishable_colors(nrImages); 
for i = 1:nrImages
    subplot(2,1,1);
    [myHist, xHist] = hist(allImages(i,:),settings.histN);
    plot(xHist, myHist, '.-','Color',myColors(i,:));   
    subplot(2,1,2);
    [myHist, xHist] = hist(normAllImages(i,:),settings.histN);
    plot(xHist, myHist, '.-','Color',myColors(i,:));
end



%% calculate medians

% Normalize if desired
if settings.useNormalizedImages
    theImages = normAllImages;
end

% median images (also for two subsets)
medianImageAll =  squeeze(median(theImages(:,:,:))); % squeeze removes singleton dim
if numel(set1)>1
    medianImage1 = squeeze(median(theImages(set1,:,:))); % squeeze removes singleton dim
else
    medianImage1 = squeeze(theImages(set1,:,:)); 
end
if numel(set2)>1
    medianImage2 = squeeze(median(theImages(set2,:,:))); % squeeze removes singleton dim
else
    medianImage2 = squeeze(theImages(set2,:,:)); 
end

% Create RGB version of those
%{
RBGmedianImageAll = repmat(medianImageAll,[1 1 3]);  
RBGmedianImage1 = repmat(medianImage1,[1 1 3]);  
RBGmedianImage2 = repmat(medianImage2,[1 1 3]);  
%}

%% show histogram of the three created imgs
figure; clf; hold on;
[myMedianHist1, xMedianHist1] = hist(medianImage1(:),settings.histN);
[myMedianHist2, xMedianHist2] = hist(medianImage2(:),settings.histN);
[myMedianHistAll, xMedianHistAll] = hist(medianImageAll(:),settings.histN);
plot(xMedianHist1,myMedianHist1,'b.-');
plot(xMedianHist2,myMedianHist2,'r.-');
plot(xMedianHist2,myMedianHistAll,'k.-');
xlabel('Distr. intensity after median multiple imgs');
ylabel('Count');
legend({'Uneven subset','Even subset','All imgs'},'Location','NorthWest')
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);


%% show images
figure(); clf; hold on; 
subplottight(1,3,1);
imshow(medianImage1,[]); title('based subset uneven');
subplottight(1,3,2);
imshow(medianImage2,[]); title('even subset');
subplottight(1,3,3);
imshow(medianImageAll,[]); title('median all');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);

% calculate percentage-wise-difference between subsets
deviationImage = medianImage1./medianImage2;

% create histogram of that and show user
figure; 
[myDevHist1, xDevHist1] = hist(deviationImage(:),settings.histN);
plot(xDevHist1, myDevHist1,'.-');
xlim([0.8,1.2]);
xlabel('distribution (matrix subset1)./(matrix subset2) values');
ylabel('count');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);

%% Final output   
outputImage = medianImageAll;

end
    
    





