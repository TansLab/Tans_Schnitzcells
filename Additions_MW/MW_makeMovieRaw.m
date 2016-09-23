function MW_makeMovieRaw(p)
% function MW_makeMovieRaw(p)


%% 
outputDir = [p.analysisDir 'movies\raw\'];
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% Loop over currently set frames
for frameIndex = settings.currentFrameRange
    
    % Load segmentation from this frame
    load([p.segmentationDir 'pos2cropseg' sprintf('%03d', frameIndex) '.mat'],'Lc','rect');
    
    % If checked segmentation exists
    if exist('Lc','var')      
        
        % print it using PN_imshowlabel
        h2=figure(); clf;
        PN_imshowlabel(p,Lc,rect,[],[]);

        % Save that figure
        saveas(h2, [outputDir 'frame' sprintf('%03d', frameIndex) '.tif']);
    else
        
        % If no seg found, skip this frame
        disp(['Frame ' num2str(frameIndex) ' was not checked, skipped for movie.']);
    end
   
    clear Lc;
    
end

end