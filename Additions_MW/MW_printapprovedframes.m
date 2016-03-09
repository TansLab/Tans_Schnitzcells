
function [approved_frame_nrs, all_frame_nrs] = MW_printapprovedframes(p)
% Print whitelist of files.
% MW 2014/09/01
%
% Loops over segmentation files in <p.segmentationDir>, loads the files,
% and checks whether Lc is precent. If so, adds that framenr to the list.
% If Lc is present, this indicates that that the segmentation in that 
% frame has been approved. Thus, a whitelist of frames is created.
%
% Input: p
%
% Output:
% - approved_frame_nrs: all approved frame nrs.
% - all_frame_nrs: all numbers of the segmentation file
%
% function [approved_frame_nrs, all_nrs] = MW_printapprovedframes(p)

% CONFIGURATION SETTINGS
FRAMENRLENGTH = 3; % length of numbers in naming. (e.g. 4 -> 004 if FRAMENRLENGTH = 3)

% Get list of segmentation files (same as .._manualcheckseg).
% ===

% Set name movie
if ~existfield(p,'outprefix')
    p.outprefix = [p.movieName 'seg'];
end

% Get files, sort
D = dir([p.segmentationDir, p.outprefix, '*.mat']);
[s,I] = sort({D.name}');
D = D(I);

% Determine position in name where frame nr is given
numposStart= length(p.outprefix)+1;
numposEnd  = findstr(D(1).name, '.mat')-1;

% Get frame nrs
segNameStrings = char(s);
all_frame_nrs = str2num(segNameStrings(:,numposStart:numposEnd))';

% Loop over segmenation files by number and check whether Lc is in there
approved_frame_nrs = [];
for i = all_frame_nrs
    clear Lc;
    
    % load file for frame #i
    load([p.segmentationDir, p.outprefix, sprintf(['%0',num2str(FRAMENRLENGTH),'d'], i), '.mat']);
    
    % If Lc is in there, frame was approved, and can be added to list of
    % approved frames
    if exist('Lc')
        approved_frame_nrs = [approved_frame_nrs i];
    end
end

end


