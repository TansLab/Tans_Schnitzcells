

function p = MW_calculateSchnitzPropertiesWithoutTracking(p, varargin)
% This function is a copy of MW_tracker, but does not perform tracking.
% This is useful because sometimes you want to retrack a few frames using
% p.overwrite=1 and MW_tracker (or another tracker). If you do this,
% properties that are required to create the schnitz are calculated only
% for the frames of interest.
% Hence, to update the schnitzFile, also the rest needs to be updated. 
% (Note: of course, the structure of these functions could be edited a bit 
% so this process would go more smoothly. But this is a TODO.)
% In any case, this function allows you to loop over all frames once more,
% and calculate the desired properties such that you can create a
% .mat schnitzcells file based on all frames.
% (This function is called from the masterscript.)
%
% Again a note: this can obviously all be done more elegantly, as this
% whole function is a direct copy+paste from MW_tracker, but with the
% actual tracking teared out.
% ***************************************************************************
%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
%% When executed as function:
numRequiredArgs = 1;
functionName = 'MW_tracker';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
  error([errorMessage]);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end

%% Processing
% lineageName is the primary tracking output
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

% by default we track only checked frames, but by setting this flag 
% you can track the uncorrected segmentation 
if ~existfield(p,'trackUnCheckedFrames')
    p.trackUnCheckedFrames = 0;
end

if ~existfield(p,'overwrite')
  p.overwrite = 0;
end
if existfield(p,'override') % backwards compatibility
  p.overwrite = p.override;
  disp('Please use p.overwrite instead of p.override');
end

% Get names of segmentation files in segmentation directory
outprefix = [p.movieName 'seg'];
D = dir([p.segmentationDir, outprefix, '*.mat']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;

% If explicit manualRange is not given, take all segmentation files
if ~existfield(p,'manualRange')
  segNameStrings = char(S);
  p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
end

% Keep only the frames in the range that contain a corrected segmentation (unless we're tracking un-checked frames)
manualRangeChecked = [];
for frameNum = p.manualRange
  clear Lc 
  load([p.segmentationDir,p.movieName,'seg',str3(frameNum)]);
  if exist('Lc','var')==1 | p.trackUnCheckedFrames 
    manualRangeChecked = [manualRangeChecked frameNum];
  else
    disp(['Skipping frame ' str3(frameNum) ' (segmentation not corrected). Use p.trackUnCheckedFrames=1 to track unchecked frames.']);
  end
end
p.manualRange = manualRangeChecked;

% if no frames in manualRange, exit here
if length(p.manualRange)==0
  error('No frames found to track. Use p.trackUnCheckedFrames=1 to track unchecked frames.');
end

disp(['Tracking ' num2str(length(p.manualRange)) ' frames ', num2str(p.manualRange(1)), ' to ', num2str(p.manualRange(end))]);
%--------------------------------------------------------------------------


%% Martijn's tracking
count = 2;  % because core schnitz functions ignore frame labels, 
            % and always start counting at 1.
for i = 2:numel(p.manualRange)
    
    % frame index
    frameIdx = p.manualRange(i);
   
    % Actual tracking MW --------------------------------------------------
    disp(['Inspecting frame ' num2str(p.manualRange(i)) '.']);
    myFileStringStart = [p.dateDir p.movieName '\segmentation\' p.movieName 'seg'];
    segFile1Path = [myFileStringStart  sprintf('%03d', p.manualRange(i-1)) '.mat'];
    segFile2Path = [myFileStringStart  sprintf('%03d', p.manualRange(i)) '.mat'];
    % End actual tracking MW ----------------------------------------------
            
    % Some stats required for checking the tracking later:
    % First frame only
    if count==2
        [Lc_fullsize_centered, Lc_fullsize, Lc] = MW_loadLcData(segFile1Path); % TODO redundancy with my fn above!
        rp  = regionprops(Lc_fullsize,'Centroid','Orientation','MajorAxisLength','Area');
        rp2 = regionprops(Lc_fullsize_centered,'Centroid');
        num_pts = size(rp,1);
        for j=1:num_pts
            Points(j).cenx      = rp(j).Centroid(1);
            Points(j).ceny      = rp(j).Centroid(2);
            Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
            Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
            Points(j).ang       = rp(j).Orientation;
            Points(j).len       = rp(j).MajorAxisLength;
            Points(j).areapx    = rp(j).Area;  %NW 2013-12
            Points(j).cellno    = j;
            Points(j).frextra   = p.manualRange(i-1); % MW DEBUG REMOVE
        end
        opts{1}=Points(1:num_pts);
    end    

    % Other frames
    % (Code from NW_tracker_centroid_vs_area)    
    [Lc_fullsize_centered, Lc_fullsize, Lc] = MW_loadLcData(segFile2Path); % TODO redundancy with my fn above!
    rp  = regionprops(Lc_fullsize,'Centroid','Orientation','MajorAxisLength','Area');
    rp2 = regionprops(Lc_fullsize_centered,'Centroid');
    num_pts = size(rp,1);
    for j=1:num_pts
        Points(j).cenx      = rp(j).Centroid(1);
        Points(j).ceny      = rp(j).Centroid(2);
        Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
        Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
        Points(j).ang       = rp(j).Orientation;
        Points(j).len       = rp(j).MajorAxisLength;
        Points(j).areapx    = rp(j).Area;  %NW 2013-12
        Points(j).cellno    = j;
        Points(j).frextra   = p.manualRange(i); % MW DEBUG REMOVE
    end
    opts{count}=Points(1:num_pts);   
    
    count = count+1;
       
    % Skip this pair if tracking file is newer than segfile
    % (This is decided in MW_linkframes, with linklistschnitz==0 as flag.)
    %{
    % (opts do need to be recalculated, so shouldn't do skipping.)
    if linklistschnitz==0, continue, end;
    %}
    
end

%% Save to schnitz format (in .mat)
% Reculculating whole lineage file from tracking files..
% (No idea how following code works, 
% stolen from NW_tracker_centroid_vs_area)
% -MW

if count>2 % if frames tracked at all
    
    % Reculculating whole lineage file from tracking files..
    MW_makeSchnitzFileFromTracking(p, opts);
    
else
    disp('WARNING: Didn''t track anything! Correct segfiles, set p.overwrite=1 to track frames if desired.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


