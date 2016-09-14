

% Extension written by Martijn Wehrens.
% Load the appropriate .mat file with skeleton data:
% data/pos1crop-skeletonData.mat 

%% Load segmentation files
load ([p.segmentationDir p.movieName 'seg' sprintf('%03d',framenr) '.mat']);

%%
theOutputdir = [p.analysisDir 'straightenedCells\skeletons\'];

if ~exist(theOutputdir,'dir')
    mkdir(theOutputdir)
end;

%%
for i = 1:length(allExtendedSkeletons)
    for j = 1:length(allExtendedSkeletons{i})

    h2=figure(); clf;
    axis equal; hold on;
    
    plot(allExtendedSkeletons{i}{j}(:,1),allExtendedSkeletons{i}{j}(:,2));    
    plot(allExtendedSkeletons{i}{j}(:,1),allExtendedSkeletons{i}{j}(:,2));
    plot(allEdges{i}{j}(:,1),allEdges{i}{j}(:,2));
    
    saveas(h2,[theOutputdir 'fr' num2str(i) 'cellno' num2str(j) '.tif']);
    
    end
end

