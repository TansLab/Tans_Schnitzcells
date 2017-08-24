function outputImage = MW_background_averaging(myFolder, myImages, ourSettings)
% function outputImage = MW_background_averaging(myFolder, myImages, ourSettings)
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
% - ourSettings            Struct for optional ourSettings. Default values
%                       hard-coded.
%   ourSettings.histN                  Number of bins used for histogram
%   ourSettings.maxImages              Matlab doesn't handle loading a large nr of images
%                                   very well. So this is capped!
%   ourSettings.neighborAveraging      Sets degree of smoothing. (How many neighbors to
%                                   take into account for neighborhood averaging of 
%                                   pixels.) Only use for shading, not for
%                                   flatfield!
%   ourSettings.useNormalizedImages    Whether imgs should be normalized.
%   ourSettings.backgroundImage        If a (same sized) background image
%                                   (sometimes called 'flatfield' throughout the script)
%                                   is given, this will be substracted from
%                                   the input images.
%
% Example list of input parameters and calling function:
% %%% Background image
%   myFolder = 'F:\A_Tans0_step1_incoming_not_backed_up\2015-05-27_background\OPT-flatfield\10ms\'
%   myImages = 'auto'; 
%   ourSettings.neighborAveraging = 0;
%   ourSettings.useNormalizedImages = 0;
%   ourSettings.maxImages = 500;
%   ourSettings.histN = 50;
%   outputImage = MW_background_averaging(myFolder, myImages, ourSettings);
%
%
% Output is img that can be used as shading image, shuold be saves as
% .mat file.
%%%%%%%%%%%%%%%%%%

% Parameter settings
if ~isfield(ourSettings, 'useNormalizedImages')
    ourSettings.useNormalizedImages = 1; % whether normalized imgs should be used
        % MW TODO OPTIMALIZATION (optional) note that code creates 
        % normalized version nevertheless
    warning('useNormalizedImages set to default value 1');
end
if ~isfield(ourSettings, 'histN')
    ourSettings.histN = 50; % number of bins used for histogram.
    warning('histN set to default value 50');
end
if ~isfield(ourSettings, 'maxImages')
    ourSettings.maxImages = 200; % max. number of images that are used
    warning('maxImages set to default value 200');    
end
if ~isfield(ourSettings, 'neighborAveraging')
    ourSettings.neighborAveraging = 5; % 
    warning('neighborAveraging set to default value 5');    
end 


%% Loading and processing images

% Automatically create list images if desired
if ~iscell(myImages), if (myImages == 'auto')
    myImages = dir(myFolder);     
    myImages = {myImages(3:end).name}; % get names and filter out '.' and '..'
end, end

% cap number of images used
if numel(myImages) > ourSettings.maxImages
    myImages = myImages(1:ourSettings.maxImages); 
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
    if isfield(ourSettings, 'backgroundImage')
        currentImage = currentImage-ourSettings.backgroundImage;
        disp('Substracting background.');
    end
    % store raw image
    allImages(i,:,:) = currentImage;
    % store filtered and normalized image    
    currentImage = currentImage/median(currentImage(:)); % normalizing by median works better than max    
    if ourSettings.neighborAveraging>0 % 
        % smooth
        currentImage = medfilt2(currentImage,[ourSettings.neighborAveraging,ourSettings.neighborAveraging]); 
    end
    normAllImages(i,:,:) = currentImage;
    % update user
    disp(['Loading ' num2str(i) '/' num2str(nrImages) '.']);
end

%% plot all histograms
h1=figure; clf; 
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
    [myHist, xHist] = hist(allImages(i,:),ourSettings.histN);
    plot(xHist, myHist, '.-','Color',myColors(i,:));   
    subplot(2,1,2);
    [myHist, xHist] = hist(normAllImages(i,:),ourSettings.histN);
    plot(xHist, myHist, '.-','Color',myColors(i,:));
end



%% calculate medians

% Use normalized images if desired
if ourSettings.useNormalizedImages
    theImages = normAllImages;
else
    theImages = allImages;
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
h2=figure; clf; hold on;
[myMedianHist1, xMedianHist1] = hist(medianImage1(:),ourSettings.histN);
[myMedianHist2, xMedianHist2] = hist(medianImage2(:),ourSettings.histN);
[myMedianHistAll, xMedianHistAll] = hist(medianImageAll(:),ourSettings.histN);
plot(xMedianHist1,myMedianHist1,'b.-');
plot(xMedianHist2,myMedianHist2,'r.-');
plot(xMedianHistAll,myMedianHistAll,'k.-');
xlabel('Distr. intensity after median multiple imgs');
ylabel('Count');
legend({'Uneven subset','Even subset','All imgs'},'Location','NorthWest')
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);


%% show images
h3=figure(); clf; hold on; 
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
h4=figure; 
[myDevHist1, xDevHist1] = hist(deviationImage(:),ourSettings.histN);
plot(xDevHist1, myDevHist1,'.-');
xlim([0.8,1.2]);
xlabel('distribution (matrix subset1)./(matrix subset2) values');
ylabel('count');
%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal'); set(gca,'FontSize',15);

%% Final output   
outputImage = medianImageAll;

if isfield(ourSettings,'figDir') & isfield(ourSettings,'figName')
    saveas(h1,[ourSettings.figDir ourSettings.figName 'fig1.fig']); 
    saveas(h2,[ourSettings.figDir ourSettings.figName 'fig2.fig']); 
    saveas(h3,[ourSettings.figDir ourSettings.figName 'fig3.fig']); 
    saveas(h4,[ourSettings.figDir ourSettings.figName 'fig4.fig']); 
    
    saveas(h1,[ourSettings.figDir ourSettings.figName 'fig1.svg']); 
    saveas(h2,[ourSettings.figDir ourSettings.figName 'fig2.svg']); 
    saveas(h3,[ourSettings.figDir ourSettings.figName 'fig3.svg']); 
    saveas(h4,[ourSettings.figDir ourSettings.figName 'fig4.svg']);
end

end
    
    





