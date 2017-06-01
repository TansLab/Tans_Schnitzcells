


%% Script to morphologically close all cells in all frames..


%% 
disp('Be sure to make a backup of your seg files!');
pause(5);

%%
% Loop over frames
lastFrameIdx=max(ourSettings.currentFrameRange);
for frIdx = ourSettings.currentFrameRange
    
    %% 
    
    % Load segfile
    clear Lc phsub LNsub rect timestamp phaseFullSize tempsegcorrect;
    load([p.segmentationDir p.movieName 'seg' num2str(frIdx) '.mat']);    
    
    % perform the closing
    Lc=NW_imclose_eachArea(Lc,strel('disk',10));
    
    % save
    save([p.segmentationDir p.movieName 'seg' num2str(frIdx) '.mat'],'Lc','phsub','LNsub','rect','timestamp','phaseFullSize','tempsegcorrect');
    
    % Update user progress
    disp(['Closed and saved all cells in frame ' num2str(frIdx) ' (of total ' num2str(lastFrameIdx) ')']);
    
    %{
    if mod(frIdx,100)==0
        disp(['Now at frame ' num2str(frIdx) ' of ' num2str(max(ourSettings.currentFrameRange))]);
    end
    %}
    
    % figure; imshow(Lc,[]);
    
    
end
    
disp('All closing done and saved.');   

