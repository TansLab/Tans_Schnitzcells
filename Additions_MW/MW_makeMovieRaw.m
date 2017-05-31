function MW_makeMovieRaw(p,ourSettings)
% function MW_makeMovieRaw(p)
%
% Note that there are some extra display option
% 
% E.g.:
% p.showPerim = 1       display outlines of cells
% p.showNr = 0,1,2      display cell numbers frame (=1) or schnitz nrs (=2)%                       
% p.slookup             give lookup table to display schnitz nrs.
% p.showPhaseImage      display phase image
% p.customColors        can be used to display schnitz colors
% p.problemCells        highlights cells in checkerboard pattern

%% 
outputDir = [p.analysisDir 'movies\raw\'];
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% add a color map if schnitz number lookup table is available
if isfield(p,'slookup')
    p=MW_addschnitzcolorscustomcolormap(p);
end

%% Loop over currently set frames
for frameIndex = ourSettings.currentFrameRange
    
    p.currentFrame = frameIndex;
    
    % Load segmentation from this frame
    load([p.segmentationDir p.movieName 'seg' sprintf('%03d', frameIndex) '.mat']);%,'Lc','rect');
    
    % If checked segmentation exists
    if exist('Lc','var')      
        
        % print it using PN_imshowlabel
        h2=figure(2); %clf; already done in PN_imshowlabel if figure is the same
        PN_imshowlabel(p,Lc,rect,[],[],'phaseImage',phsub);

        % Save that figure
        h2.InvertHardcopy = 'off'; % some weird setting preventing post-processing that changes text color (see http://nl.mathworks.com/help/matlab/creating_plots/save-figure-preserving-background-color.html)
        saveas(h2, [outputDir 'frame' sprintf('%03d', frameIndex) '.tif'],'tif');
    else
        
        % If no seg found, skip this frame
        disp(['Frame ' num2str(frameIndex) ' was not checked, skipped for movie.']);
    end
   
    clear Lc;
    
end

end