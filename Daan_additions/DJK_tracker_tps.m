function [p, processed] = DJK_tracker_tps(p,varargin)
%   
% Same as tracker_tps with some minor changes. Added feature that if
% segmentation files are younger than last tracking between frames, the
% tracking is done again (and not skipped).
%
%   TRACKER_TPS allows users to track cells identified by cell segmentation 
%   across all frames for a movie.  Cells in one frame are matched to cells 
%   in adjacent frames, and cell division is recognized by associating each 
%   parent cell with its children.
%   
%   TRACKER_TPS(P,'Field1',Value1,'Field2',Value2,...) also performs cell
%   tracking using segmentation results, but the extra arguments permit users 
%   to adjust any parameters describing the movie or parameters controlling 
%   the cell tracking process.  The extra arguments can overwrite any specific 
%   parameters provided in P by setting P.Field1 = Value1, P.Field2 = Value2, 
%   etc.  Thus any/all schnitzcells parameter values can be defined in the 
%   function call via these optional field/value pairs.  (This is in the style 
%   of setting MATLAB properties using optional property/value pairs.)  
%   
%   TRACKER_TPS produces an output file (a MATLAB binary file) that 
%   contains the 'schnitzcells' variable describing each cell's lineage. The 
%   path name of this file defaults to [p.tracksDir,p.movieName,'_lin.mat'], 
%   but can be specified by providing an optional 'lineageName' field in the 
%   schnitcells parameter structure, or by providing an extra 
%   'lineageName',value pair of arguments to this function. 
%   
%   TRACKER_TPS returns a struct (1x1 struct array) referred to as the 
%   schnitzcells parameter structure that contains fields and values 
%   contained in P, including unchanged/original parameters plus any of those 
%   added or overridden in the list of properties & values provided in the 
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control tracking
% 
%   lineageName   the schnitzcells cell lineage structure is written to this 
%                 file, by default named [p.trackDir p.moveiName '_lin.mat']
%   
%   trackRange    range of frame numbers to track across; by default all 
%                 frames containing segmentation files with corrected 
%                 segmentation fields will be used in tracking
%   
%   trackUnCheckedFrames  flag that when set to true or 1 permits this 
%                 routine to perform tracking on segmentation files that were 
%                 not manually verified or corrected (i.e. they do not have Lc)
%   
%   overwrite      flag that when set to true or 1 permits this routine to 
%                 perform tracking on pairs of frames even if tracking was 
%                 previously performed for those frames
%   
%-------------------------------------------------------------------------------
%
% Varibles contained in the output file:
%
%  schnitzcells   schnitzcells is a 1xNS struct array (NS = number of "schnitz" 
%                 tracks).  This struct array describes the lineage of all 
%                 cells in a movie.  A "schnitz" is a track that begins at 
%                 the first frame a distint cell appears in, and ends at the 
%                 last frame that distinct cell appears in.  When a cell 
%                 divides it's track ends, and two new tracks being at the 
%                 first frame it's two children appear in.  A schnitz (track) 
%                 is numbered by its position (index) within the struct array. 
%                 A minimal schnitz structure is described by P (the schnitz 
%                 number of it's parent), children schnitz numbers D and E,
%                 sister schnitz number S, frames (an array) listing each frame 
%                 this schnitz occurs in, cellno (an array) listing each 
%                 cell number of the tracked cell within each of the frames, 
%                 and N, the number of frames that distinct cell appears.
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
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end

if ~existfield(p,'overwrite')
  p.overwrite = 0;
end
if existfield(p,'override') % backwards compatibility
  p.overwrite = p.override;
  disp('Please p.overwrite use instead of p.override');
end

% lineageName is the primary tracking output
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
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

% by default, will use partial results
%if ~existfield(p,'usepartial')
%    p.usepartial = true
%end

firstframename = [p.segmentationDir, outprefix, str3(p.trackRange(1)), '.mat'];
segData = load(firstframename);

% added nitzan 2005June24
if ~isfield(segData,'Lc') & p.trackUnCheckedFrames
    segData.Lc = segData.LNsub;
end
% end addition nitzan 2005June24
pts =     regionprops(segData.Lc,'Centroid');
orn =     regionprops(segData.Lc,'Orientation');
lngth =   regionprops(segData.Lc,'MajorAxisLength');
% Find all invalid points
Invalid = find([lngth.MajorAxisLength]==0);
if (size(Invalid,2)~=0)
    pts([Invalid]).Centroid = [0 0];
    orn([Invalid]).Orientation = 0;
    lngth([Invalid]).MajorAxisLength = 1;
end;

num_pts = size(pts,1);

for j=1:num_pts
    Points(j).cenx = pts(j).Centroid(1);
    Points(j).ceny = pts(j).Centroid(2);
    Points(j).ang = orn(j).Orientation;
    Points(j).len = lngth(j).MajorAxisLength;
    Points(j).cellno=j;
end

opts{1}=Points(1:num_pts);
  

% open up output tracks info
% seq_out_name = [p.tracksDir, p.movieName, '_tracks']; % DJK 081022
% fid_out = fopen([seq_out_name '.txt'],'w'); % DJK 081022


% before tracking pairs, get info on TPS execuable and directory 
trackerFullPath = which('tracker_tps');
if isempty(trackerFullPath)
  error('TPS tracker path not found!')
end
trackCodeDir = [fileparts(trackerFullPath) filesep];
str = computer();
switch(str)
    case 'PCWIN'
        trackExecutable = 'DotpsNW';
    case 'MAC'
        trackExecutable = 'dotps.mac';
    case 'GLNX86'
        trackExecutable = 'dotps.linux';
    otherwise
        error (['There is no tracker executable for a ' str ' computer!']);
end

% Work has to be done in the track code directory so executable 
% can find the softconfig.dat file (lame, I know)
originalDir = pwd;
cd(trackCodeDir)

% processed will contain numbers of frames that were processed
processed = [];

% CARRY OUT TRACKING
i = 2;  % i = row index into cellmat, NOT frameNum
for frameNum = p.trackRange(2:end)
    todaySegFile= [p.segmentationDir, outprefix, str3(frameNum), '.mat'];
    mynum_t= str3(frameNum);
    yesterdayFrameNum = p.trackRange(find(p.trackRange==frameNum)-1);
    yesterdaySegFile = [p.segmentationDir, outprefix, str3(yesterdayFrameNum), '.mat']; %DJK 081107 added '.mat'
    mynum_y= str3(yesterdayFrameNum);

%    if (p.usepartial)
%        if (frameNum >= lastFrameNumTracked)
%            disp(['Skipping previously tracked frame pair ',...
%		    mynum_y,' , ',mynum_t]);
%	    i = i+1;
%	    continue
%        end
%    end

    fprintf(['Starting to process frame pair ',mynum_y,' , ',mynum_t]);

    yesterdayTrackInputFile = [p.tracksDir,p.movieName,'-tps-input-',...
                               mynum_y,'.txt'];
    if (i==2)
        % only write yesterday's data the first time processing the pair
        % why? because today becomes tomorrow's yesterday, right?
        yesterdaySegData = load(yesterdaySegFile);
        % added nitzan 2005June24
        if ~isfield(yesterdaySegData,'Lc') & p.trackUnCheckedFrames
            yesterdaySegData.Lc = yesterdaySegData.LNsub;
        end
        % end addition nitzan 2005June24
        
        write_tps_tracker_inputs(yesterdaySegData.Lc,yesterdayTrackInputFile);
    end

    todayTrackInputFile = [p.tracksDir,p.movieName,'-tps-input-',...
                           mynum_t,'.txt'];
    todaySegData = load(todaySegFile);
    % added nitzan 2005June24
    if ~isfield(todaySegData,'Lc') & p.trackUnCheckedFrames
        todaySegData.Lc = todaySegData.LNsub;
    end
    % end addition nitzan 2005June24

    write_tps_tracker_inputs(todaySegData.Lc,todayTrackInputFile);

    pts =     regionprops(todaySegData.Lc,'Centroid');
    orn =     regionprops(todaySegData.Lc,'Orientation');
    lngth =   regionprops(todaySegData.Lc,'MajorAxisLength');
    % Find all invalid points
    Invalid = find([lngth.MajorAxisLength]==0);
    if (size(Invalid,2)~=0)
        for djk = Invalid %DJK 071207
            pts(djk).Centroid = [0 0]; %DJK 071207
            orn(djk).Orientation = 0; %DJK 071207
            lngth(djk).MajorAxisLength = 1; %DJK 071207
        end %DJK 071207
    end;
    
    num_pts = size(pts,1);

    for j=1:num_pts
        Points(j).cenx = pts(j).Centroid(1);
        Points(j).ceny = pts(j).Centroid(2);
        Points(j).ang = orn(j).Orientation;
        Points(j).len = lngth(j).MajorAxisLength;
        Points(j).cellno=j;
    end

    opts{i}=Points(1:num_pts);

    logFile         = [p.tracksDir,p.movieName,'-tps-log-',...
                       mynum_y,'-to-',mynum_t,'.txt'];

    trackOutputFile = [p.tracksDir,p.movieName,'-tps-output-',...
                       mynum_y,'-to-',mynum_t,'.txt'];

    trackCommand    = [trackCodeDir,trackExecutable,' ',...
                       '"', yesterdayTrackInputFile,'" "',...
                       todayTrackInputFile,'"']; % DJK 071129
               
    % Redirect the output of the tracking command to a log file
    str = computer();
    switch (str)
        case 'PCWIN'
            trackCommand = [trackCommand,' > "',logFile,'" 2>&1']; % DJK 071129
        otherwise
            trackCommand = [trackCommand,' >& "',logFile,'"']; % DJK 071129
    end
 
% START DJK 081107
% if trackOutputFile is younger than segmentation files, means that
% segmentation has been updated, and we want to track again
    segUpdated = false;
    if exist(trackOutputFile)==2
      info_yesterdaySegFile = dir(yesterdaySegFile);
      info_todaySegFile = dir(todaySegFile);
      info_trackOutputFile = dir(trackOutputFile);
      datenumber_yesterdaySegFile = datenum(info_yesterdaySegFile.date);
      datenumber_todaySegFile = datenum(info_todaySegFile.date);
      datenumber_trackOutputFile = datenum(info_trackOutputFile.date);
      if datenumber_yesterdaySegFile>datenumber_trackOutputFile | datenumber_todaySegFile>datenumber_trackOutputFile
        segUpdated = true;
      end
    end
% END DJK 081107
   
    if exist(trackOutputFile)==2 & ~(p.overwrite) & ~segUpdated% DJK 081107
      fprintf(1,' -> Skipping cause previously tracked (use p.overwrite=1 to redo)\n');
      i = i+1;
      continue
    end
    fprintf(1,'\n');

    % clean up any existing output
    if exist('Match.dat') == 2
        delete('Match.dat')
    end
    if exist(trackOutputFile) == 2
        delete(trackOutputFile)
    end

    % execute track command
    system(trackCommand);

    % move cell matching results to track output file
    [status,msg] = movefile('Match.dat',trackOutputFile);
    if (status == 0)
        disp(['Warning: encountered problem matching frame pair ',...
              mynum_y,' , ',mynum_t,'...']);
    end
    
    i = i+1;
    
    % update processed frames
    processed = [processed yesterdayFrameNum frameNum];
end
    
% get unique processed frames
processed = unique(processed);
    
% Convert matching results to schnitzcells-format lineage
[P D E G]= data_treat(p);


if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin_tps.mat'];
end

[schnitzcells cellid]=recalc_schnitz(P,D,E,G,p.trackRange,'',opts);
trackRange = p.trackRange;

schnitzcells = renumberschnitzes(p,schnitzcells);

disp(['saving schnitzcells lineage structure to ' p.lineageName]);
save(p.lineageName,'schnitzcells');

% Go back to original directory now that tracking is complete
cd(originalDir)

