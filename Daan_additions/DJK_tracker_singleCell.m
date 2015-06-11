function [p] = DJK_tracker_singleCell (p, varargin)
% DJK_SINGLECELLTRACK   Make a TrackComplete output, but then without tracking.
% 
%   Copied and Adjusted from TRACKCOMPLETE 
%   Has same beginning as TRACKCOPLETE
%   
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------
numRequiredArgs = 1;
if (nargin < 1) | ...
   (mod(nargin,2) == 0) | ...
   (~isSchnitzParamStruct(p))
  errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> trackcomplete:',...
      '    Invalid input arguments.',...
      '    Try "help trackcomplete".');
  error(errorMessage);
end

%-------------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values 
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
          'Error using ==> trackcomplete:',...
          '    Invalid property ', num2str(varargin{i}), ...
          ' is not (needs to be) a string.',...
          '    Try "help trackcomplete".');
      error(errorMessage);
    end
    fieldName = schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end

% lineageName is the primary tracking output
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

% trackName file holds the original tracker's cellmat output
if ~existfield(p,'trackName')
  p.trackName = [p.tracksDir,p.movieName,'_Tdata.mat'];
end

% by default we track only checked frames, but by setting this flag 
% you can track the uncorrected segmentation 
if ~existfield(p,'trackUnCheckedFrames')
    p.trackUnCheckedFrames = 0;
end

% Get directory of existing segmentation files 
outprefix = [p.movieName 'seg'];
D = dir([p.segmentationDir, outprefix, '*.mat']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;

% If explicit track range not given, figure it out by determining 
% which segmented frames have been corrected
if ~existfield(p,'trackRange')
  D = dir([p.segmentationDir, outprefix, '*.mat']);
  [S,I] = sort({D.name}');
  D = D(I);
  frameNumStrings = char(S);
  p.trackRange = str2num(frameNumStrings(:,numpos:numpos+2))';
end

% Keep only the frames in the track range that contain a corrected segmentation
% (unless we're tracking un-checked frames)
segmentedFrames = [];
for frameNum=p.trackRange
  clear Lc;
  tempsegcorrect=0;
  load([p.segmentationDir,p.movieName,'seg',str3(frameNum)]);
  if (exist('Lc')==1 & tempsegcorrect==0) | p.trackUnCheckedFrames 
    segmentedFrames = [segmentedFrames frameNum];
  else
    disp(['Skipping frame ' num2str(frameNum) ... 
	    ' (segmentation not corrected). Use p.trackUnCheckedFrames=1 to track unchecked frames.']);
  end
end
p.trackRange=segmentedFrames;

if length(p.trackRange)==0
    error('No frames found to track. Use p.trackUnCheckedFrames=1 to track unchecked frames.');
end
lastname= [p.segmentationDir, outprefix, str3(p.trackRange(end)), '.mat'];
clear LNsub Lc;
load(lastname);
if exist('Lc') 
    MaxCell= max2(Lc);
elseif p.trackUnCheckedFrames   % added nitzan 2005June24
    MaxCell= max2(LNsub);       % added nitzan 2005June24
else
    error(['segmentation for frame ' num2str(p.trackRange(end)) ...
	    ' not verified/corrected.']);
    MaxCell= max2(LNsub);
end;
tcells= 1:MaxCell;
cellmat= zeros(length(p.trackRange), MaxCell);
cellmat(length(p.trackRange),:)= 1:MaxCell;
se= strel('diamond', 1);






% -------------------------------------------------------------------------
% From here on: DJK 071023

i = 1;
for frameNum = 1:length(p.trackRange) % Loop over frames in segRange

    % get the segmentation data from the current frame
    framename = [p.segmentationDir, outprefix, str3(p.trackRange(frameNum)), '.mat'];
    segData = load(framename);
    if ~isfield(segData,'Lc') & p.trackUnCheckedFrames
        segData.Lc = segData.LNsub;
    end
    pts =     regionprops(segData.Lc,'Centroid');
    orn =     regionprops(segData.Lc,'Orientation');
    lngth =   regionprops(segData.Lc,'MajorAxisLength');
    
    disp(['Getting data from ' framename]);

    for cell = 1:length(pts) % Loop over de the cells in the current frame
        % Get data from current cell and save as part of temp
        temp.P = 0;
        temp.E = 0;
        temp.D = 0;
        temp.frame_nrs = p.trackRange(frameNum); %DJK 080704 changed from frameNum+1; // MW fix N+1 bug 2014/06/24
        temp.cenx = pts(cell).Centroid(1);
        temp.ceny = pts(cell).Centroid(2);
        temp.ang = orn(cell).Orientation;
        temp.len = lngth(cell).MajorAxisLength;
        temp.cellno = cell;
        test(i) = temp;
        i = i+1;
    end;
end;

% Save temp as an object called "schnitzcells" to a file
schnitzcells = test;
disp(['saving schnitzcells lineage structure to ' p.lineageName]);
save(p.lineageName,'schnitzcells');
