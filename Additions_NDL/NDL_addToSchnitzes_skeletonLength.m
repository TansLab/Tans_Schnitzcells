function [p, schnitzcells] = NDL_addToSchnitzes_skeletonLength(p)
% load schnitzcells
load([p.tracksDir p.movieName '-Schnitz.mat'])

% calculate all lengths
[p, allLengthsOfBacteriaInPixels, allLengthsOfBacteriaInMicrons] = NDL_lengthforfillamentedcells(p, [1:10])

%assign to schnitz
%(loops)
% use schnitzcells.frame_nrs and schnitzcells.cellno

%% save
save([p.tracksDir p.movieName '-Schnitz.mat'])
end