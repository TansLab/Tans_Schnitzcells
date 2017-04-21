

% In the newer code, the field 'frames' is replaced by the field
% 'frame_nrs', because frame_nrs has the correct numbering, and 'frames'
% has a value that is 1 too high. See also MW_calculateframe_nrs.



load([p.tracksDir p.movieName '-Schnitz.mat']);

schnitzcells = MW_calculateframe_nrs(schnitzcells)

save([p.tracksDir p.movieName '-Schnitz.mat'], 'schnitzcells');