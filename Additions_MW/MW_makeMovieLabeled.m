function MW_makeMovieLabeled(p, ourSettings)
% function MW_makeMovieRaw(p)

% OPTIONAL ARGUMENTS
% p.showPerim=1       if p.showPerim is valid field creates perimeter movie

%% 
outputDir = [p.analysisDir 'movies\MWlabeled\'];
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% Create coloring for labels
p.slookup = MW_makeslookup(p)

highestSchnitzIndx = max(p.slookup(size(p.slookup,1),:));
% Set up custom colormap
% easy way:
%theColorMap = linspecer(maxCellNo);

% let's use one of the standard colormaps
standardColorMap = hsv(highestSchnitzIndx); % hsv jet

% but mix it up such that neighbors have different
% colors                    
shuffle=randperm(highestSchnitzIndx); 
%shuffle = [10  20  28  35  29  27  64  40  33  37  47  58  22  36  18  21  50  57  34  25  19  43   1  49  16  60  23   3  48   9  45  38  44  46   4  26   7  15  54  59  55  24   6  14  42  13   8  61  63  41   5   2  30  53  31  52  32  51  56  17  11  62  12  39];

% perform shuffling
standardColorMapShuffled = NaN(size(standardColorMap));
standardColorMapShuffled(:,1) = standardColorMap(shuffle,1);
standardColorMapShuffled(:,2) = standardColorMap(shuffle,2);
standardColorMapShuffled(:,3) = standardColorMap(shuffle,3);

% how many copies of the colormap do we need?
%copiesNeeded = ceil(maxCellNo/COLORMAPSIZE);

% create the color map
%theColorMap = repmat(standardColorMapShuffled,copiesNeeded,1);

% set the colormap
p.customColors = [0 0 0; standardColorMapShuffled; 1 1 1];

%% Loop over currently set frames
for frameIndex = ourSettings.currentFrameRange
    
    %%         
    
    % Load segmentation from this frame
    load([p.segmentationDir 'pos2cropseg' sprintf('%03d', frameIndex) '.mat'],'Lc','rect','phsub');
    %load([p.segmentationDir 'pos2cropseg' sprintf('%03d', frameIndex) '.mat']);
    
    % If checked segmentation exists
    if exist('Lc','var')      
        
        % print it using PN_imshowlabel
        h2=figure(1); clf;
        p.currentFrame = frameIndex; p.showNr=2;
        PN_imshowlabel(p,Lc,rect,[],[],'phaseImage',phsub);

        % Save that figure
        saveas(h2, [outputDir 'frame' sprintf('%03d', frameIndex) '.tif']);
    else
        
        % If no seg found, skip this frame
        disp(['Frame ' num2str(frameIndex) ' was not checked, skipped for movie.']);
    end
   
    clear Lc;
    
end

end