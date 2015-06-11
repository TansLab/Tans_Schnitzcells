function problemCells = DJK_analyzeTracking(p, varargin);
% DJK_analyzeTracking analyses the complete tracking file. 
%
% *** Less convenient output file than in the newer version
% (DJK_analyzeTracking) but this older version can handle skips in frames
% (e.g. one frame was defocused and is not used). To be updated. ***

% Reports on :
% # cells tracked in complete lineage 
%   # cells in manualRange
%     * cells that are in segmentation but not in tree
%     * cells that are orphan (first frame always has orphans)
%     * cells that are barren (cell lineage ends before last frame)
%     * cells that move more than p.pixelsMoveDef pixels between frames 
%     * cells that grow more/less than p.pixelsLenDef pixels between frames 
%     * cells that change length more/less than 20 pixels after division 
%
% OPTIONAL ARGUMENTS:
%  'lineageName'    allows to use a tracking lineage file other than standard
%
%  'manualRange'    allows to analyze a subset of frames (standard: all framse)
%
%  'DJK_saveDir'    Directory where results will be saved. Defaults to
%                   "p.analysisDir 'tracking\' 'manualRange1_200\'"
%
%  'pixelsMoveDef'  threshold number for moving cell in pixels (default: 10)
%
%  'pixelsLenDef'   threshold numbers for growing cell in pixels. 
%                   default: [-4 6] = shrink 4 of grow more than 6 => weird
%                   Note: old length definition is used (fit of ellipse)
%

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1;
functionName = 'DJK_analyzeTracking';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
  error(errorMessage);
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
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If explicit lineageName is not given, use standard
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

% Load lineage file
if ~(exist(p.lineageName)==2)
  error(['Could not read tracking file ''' p.lineageName '''.']);
end
load(p.lineageName);

% If explicit manualRange is not given, take all frames that are in schnitzcells
if ~existfield(p,'manualRange')
  fr_maximum = -1; fr_minimum = 1000;
  for cell = 1:length(schnitzcells)
    if (max(schnitzcells(cell).frames) > fr_maximum) 
      fr_maximum = max(schnitzcells(cell).frames);
    end
    if (min(schnitzcells(cell).frames) < fr_minimum) 
      fr_minimum = min(schnitzcells(cell).frames);
    end
  end
  p.manualRange = [fr_minimum-1:fr_maximum-1]; % in schnitzcells frames are +1
end

% If explicit DJK_saveDir is not given, use standard
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'tracking' filesep 'manualRange' num2str(p.manualRange(1)) '_' num2str(p.manualRange(end)) filesep];
end

% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end

% If explicit pixelsMoveDef is not given, use 10 pixels
if ~existfield(p,'pixelsMoveDef')
    p.pixelsMoveDef = 10;
end

% If explicit pixelsMoveDef is not given, use [-1 4] pixels
if ~existfield(p,'pixelsLenDef')
    p.pixelsLenDef = [-4 6];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Open file to write results to
%--------------------------------------------------------------------------
fid = fopen([p.DJK_saveDir p.movieName '-tracking.txt'],'wt');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PREPARATION OF TRACKING DATA
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
dispAndWrite(fid, [' * Saving results in : ' p.DJK_saveDir]);
dispAndWrite(fid, [' * Loaded tracking file from : ' p.lineageName]);

firstFrame = min(p.manualRange);
lastFrame = max(p.manualRange);

% Determine for each cell in schnitzcells whether really in tracking
% indicated by inTracking (schnitzcells(1).inTracking=true)
% in schnitzcells frames are +1, so correct for that
range = p.manualRange + 1;
nrCellsInTracking = 0;
for cell = 1:length(schnitzcells)
  schnitzcells(cell).inTracking = false;
  if (intersect(range, schnitzcells(cell).frames))
    schnitzcells(cell).inTracking = true;
    nrCellsInTracking = nrCellsInTracking+1;
  end
end

% Get number of segmented cells of each frame
for frameNum = p.manualRange
  clear LNsub Lc timestamp;
  load([p.segmentationDir, p.movieName, 'seg', str3(frameNum)]);
  if ~exist('Lc') 
    Lc = LNsub;
  end
  NrCellsPerFrame(frameNum) = max2(Lc);
end;

% Get number of cells in tracking for each frame
NrCellsInTrackingPerFrame = zeros(lastFrame, 1);
for schnitznum = 1:length(schnitzcells)
  for frame = schnitzcells(schnitznum).frames - 1;
    if (frame<=lastFrame)
      NrCellsInTrackingPerFrame(frame) = NrCellsInTrackingPerFrame(frame) + 1;
    end;
  end
end

% Find cells that are in segmentation but not in tree
segAndTrackDoNotMatch = zeros(lastFrame, 1);
for frameNum = p.manualRange
  if ~(NrCellsPerFrame(frameNum) == NrCellsInTrackingPerFrame(frameNum))
    segAndTrackDoNotMatch(frameNum) = 1;
  end
end

% Find orphans
nrOrphansInManualRange = 0;
for cell = 1:length(schnitzcells)
  schnitzcells(cell).orphan = true;
  if schnitzcells(cell).P > 0;
    schnitzcells(cell).orphan = false;
  else
    if schnitzcells(cell).inTracking & schnitzcells(cell).frames(1) > firstFrame+1 % cells in first frame do not count
      nrOrphansInManualRange = nrOrphansInManualRange + 1;
    end
  end
end

% output of problem cells (schnitz nr + frame)
problemCells = [];

% Find barren cells
framesWithBarrenCells = [];
%nrBarrenInManualRange = 0;
for cell = 1:length(schnitzcells)
  schnitzcells(cell).barren = false;
  if schnitzcells(cell).E==0 | schnitzcells(cell).D==0;
    schnitzcells(cell).barren = true;
    if schnitzcells(cell).inTracking & schnitzcells(cell).frames(end)-1 < lastFrame % cells in last frame do not count
      framesWithBarrenCells = [framesWithBarrenCells schnitzcells(cell).frames(end)-1];
      problemCells = [problemCells ; cell schnitzcells(cell).frames(end)-1];
      % nrBarrenInManualRange = nrBarrenInManualRange + 1;
    end
  end
end
framesWithBarrenCells = unique(framesWithBarrenCells);

% Find cells moving > p.pixelsMoveDef pixels
framesWithCellsMoving = [];
for cell = 1:length(schnitzcells)
  schnitzcells(cell).moving = [];
  cenx = schnitzcells(cell).cenx_cent(1); % DJK 090410 cenx(i);
  ceny = schnitzcells(cell).ceny_cent(1); % DJK 090410 ceny(i);
  for i = 2:length(schnitzcells(cell).frames);
    if schnitzcells(cell).inTracking & schnitzcells(cell).frames(i)-1 <= lastFrame
      cenx_new = schnitzcells(cell).cenx_cent(i); % DJK 090410 cenx(i);
      ceny_new = schnitzcells(cell).ceny_cent(i); % DJK 090410 ceny(i);
      if sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 ) > p.pixelsMoveDef
        schnitzcells(cell).moving = [schnitzcells(cell).moving schnitzcells(cell).frames(i)-1];
        framesWithCellsMoving = [framesWithCellsMoving schnitzcells(cell).frames(i)-1];
        problemCells = [problemCells ; cell schnitzcells(cell).frames(i)-1];
%         disp([str3(cell) ' moved ' num2str(sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 )) ' pixels ']);
      end
      cenx = cenx_new; ceny = ceny_new;
    end
  end
end
framesWithCellsMoving = unique(framesWithCellsMoving);

% Find cells growing < > p.pixelsLenDef pixels
framesWithCellsGrowingTooLittle = [];
framesWithCellsGrowingTooMuch = [];
for cell = 1:length(schnitzcells)
  schnitzcells(cell).growingTooLittle = [];
  schnitzcells(cell).growingTooMuch = [];
  len = schnitzcells(cell).len(1);
  for i = 2:length(schnitzcells(cell).frames);
    if schnitzcells(cell).inTracking & schnitzcells(cell).frames(i)-1 <= lastFrame
      len_new = schnitzcells(cell).len(i);
      if (len_new-len) < p.pixelsLenDef(1)
        schnitzcells(cell).growingTooLittle = [schnitzcells(cell).growingTooLittle schnitzcells(cell).frames(i-1)-1];
        framesWithCellsGrowingTooLittle = [framesWithCellsGrowingTooLittle schnitzcells(cell).frames(i-1)-1];
        problemCells = [problemCells ; cell schnitzcells(cell).frames(i-1)-1];
      elseif (len_new-len) > p.pixelsLenDef(2) 
        schnitzcells(cell).growingTooMuch = [schnitzcells(cell).growingTooMuch schnitzcells(cell).frames(i-1)-1];
        framesWithCellsGrowingTooMuch = [framesWithCellsGrowingTooMuch schnitzcells(cell).frames(i-1)-1];
        problemCells = [problemCells ; cell schnitzcells(cell).frames(i-1)-1];
      end
      len = len_new;
    end
  end
end
framesWithCellsGrowingTooLittle = unique(framesWithCellsGrowingTooLittle);
framesWithCellsGrowingTooMuch = unique(framesWithCellsGrowingTooMuch);
framesWithCellsGrowingWeird = unique([framesWithCellsGrowingTooLittle framesWithCellsGrowingTooMuch]);

% Find cells that change length more/less than 20 pixels after division
framesWithCellsChangingAfterDivision = [];
for cell = 1:length(schnitzcells)
  schnitzcells(cell).offspringToBig = false;
  D = schnitzcells(cell).D;
  E = schnitzcells(cell).E;
  if schnitzcells(cell).inTracking & E>0 & D>0
    combined_length = schnitzcells(D).len(1) + schnitzcells(E).len(1);
    if (schnitzcells(cell).len(end) >  20 + schnitzcells(D).len(1) + schnitzcells(E).len(1)) | ...
       (schnitzcells(cell).len(end) < -20 + schnitzcells(D).len(1) + schnitzcells(E).len(1))
      schnitzcells(cell).offspringToBig = true;
      framesWithCellsChangingAfterDivision = [framesWithCellsChangingAfterDivision schnitzcells(cell).frames(end)-1];
      problemCells = [problemCells ; cell schnitzcells(cell).frames(end)-1];
    end
  end
end
framesWithCellsChangingAfterDivision = unique(framesWithCellsChangingAfterDivision);

% output of problem cells
problemCells = unique(problemCells, 'rows');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% CHECKING OF TRACKING
%--------------------------------------------------------------------------
% Display nr of cells in lineage file
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Results:']);
dispAndWrite(fid, ['--- ' str3(length(schnitzcells)) ' cells in tracking file']);
dispAndWrite(fid, ['  | ']);
dispAndWrite(fid, ['  |--- ' str3(nrCellsInTracking) ' of these are in manualRange']);
dispAndWrite(fid, ['     | ']);

% Display cells that are in segmentation but not in tree
dispAndWrite(fid, ['     |--- ' str3(sum(segAndTrackDoNotMatch)) ' frames lack cells in tracking, which are there in segmentation']);
for frameNum = p.manualRange
  if segAndTrackDoNotMatch(frameNum)
    dispAndWrite(fid, ['     |  |--- frame ' str3(frameNum) ' has ' str3(NrCellsInTrackingPerFrame(frameNum)) ' / ' str3(NrCellsPerFrame(frameNum)) ' cells in tracking / segmentation']);
    for cellnum = [1:NrCellsPerFrame(frameNum)]
      if ( findSchnitz(schnitzcells, frameNum+1, cellnum) == -1)
        dispAndWrite(fid, ['     |     |--- cellno ' str3(cellnum) ' is missing in tracking']);
      end
    end
  end
end

% Display orphans (first frame always has orphans)
dispAndWrite(fid, ['     |  ']);
dispAndWrite(fid, ['     |--- ' str3(nrOrphansInManualRange) ' of cells are orphans']);
for cell = 1:length(schnitzcells)
  if schnitzcells(cell).inTracking & schnitzcells(cell).orphan & schnitzcells(cell).frames(1) > firstFrame+1 % cells in first frame do not count
    dispAndWrite(fid, ['     |  |--- schnitz ' str3(cell) ' is an orphan and appears in frame ' str3(schnitzcells(cell).frames(1)-1)]);
  end
end

% Display barren cells (cell lineage ends)
dispAndWrite(fid, ['     |  ']);
% dispAndWrite(fid, ['     |--- ' str3(nrBarrenInManualRange) ' of barren cells in manualRange']);
dispAndWrite(fid, ['     |--- ' str3(length(framesWithBarrenCells)) ' frames have barren cells']);
for fr = framesWithBarrenCells
  out = ['     |  |--- in frame ' str3(fr) ' schnitz '];
  for cell = 1:length(schnitzcells)
    if schnitzcells(cell).inTracking & schnitzcells(cell).barren & intersect(fr, schnitzcells(cell).frames(end)-1)
      out = [out str3(cell) ' '];
      %dispAndWrite(fid, ['     |  |--- schnitz ' str3(cell) ' is barren and disappears in frame ' str3(schnitzcells(cell).frames(end)-1)]);
    end
  end
  out = [out 'are barren'];
  dispAndWrite(fid, out);
end

% Display moving cells (possibly wrong tracking)
dispAndWrite(fid, ['     |  ']);
dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsMoving)) ' frames have cells moving > ' num2str(p.pixelsMoveDef) ' pixels']);
for fr = framesWithCellsMoving
%  out = ['     |  |--- in frame ' str3(fr) ' schnitz '];
  dispAndWrite(fid, ['     |  |--- in frame ' str3(fr)]);
  for cell = 1:length(schnitzcells)
    if schnitzcells(cell).inTracking & intersect(fr, schnitzcells(cell).moving)
%       out = [out str3(cell) ' '];
      fr_idx = find(schnitzcells(cell).frames==(fr+1));
      cenx = schnitzcells(cell).cenx_cent(fr_idx-1);
      ceny = schnitzcells(cell).ceny_cent(fr_idx-1);
      cenx_new = schnitzcells(cell).cenx_cent(fr_idx);
      ceny_new = schnitzcells(cell).ceny_cent(fr_idx);
      distanceMoved = sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 );
      dispAndWrite(fid, ['     |  |  |--- schnitz ' str3(cell) ' : ' num2str( round(distanceMoved) ) ' pixels']);
    end
  end
%   out = [out 'are moving'];
%   dispAndWrite(fid, out);
  dispAndWrite(fid, ['     |  |' ]);
end

% Display weird growing cells (possibly wrong tracking / segmentation)
dispAndWrite(fid, ['     |  ']);
dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsGrowingWeird)) ' frames have cells growing < ' num2str(p.pixelsLenDef(1)) ' or > ' num2str(p.pixelsLenDef(2)) ' pixels']);
for fr = framesWithCellsGrowingWeird
  dispAndWrite(fid, ['     |  |--- in frame ' str3(fr)]);
  for cell = 1:length(schnitzcells)
    idx = [find(schnitzcells(cell).growingTooLittle==fr) find(schnitzcells(cell).growingTooMuch==fr)];
    if schnitzcells(cell).inTracking & idx
      fr_idx = find(schnitzcells(cell).frames==(fr+1));
      lengthIncrease = schnitzcells(cell).len(fr_idx+1) - schnitzcells(cell).len(fr_idx);
      dispAndWrite(fid, ['     |  |  |--- schnitz ' str3(cell) ' : ' num2str( round(lengthIncrease) ) ' pixels']);
    end
  end
  dispAndWrite(fid, ['     |  |' ]);
end

% Display moving whose length change after division (possibly wrong tracking)
dispAndWrite(fid, ['     |  ']);
dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsChangingAfterDivision)) ' frames have cells with length change <> 20 pixels after division']);
for fr = framesWithCellsChangingAfterDivision
  dispAndWrite(fid, ['        |--- in frame ' str3(fr)]);
  for cell = 1:length(schnitzcells)
    if schnitzcells(cell).inTracking & schnitzcells(cell).offspringToBig & fr==schnitzcells(cell).frames(end)-1
      D = schnitzcells(cell).D;
      E = schnitzcells(cell).E;
      lengthIncrease = schnitzcells(D).len(1) +  schnitzcells(E).len(1) - schnitzcells(cell).len(end);
      dispAndWrite(fid, ['        |  |--- schnitz ' str3(cell) ' : ' num2str( round(lengthIncrease) ) ' pixels']);
    end
  end
  dispAndWrite(fid, ['        |' ]);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Close file to write results to
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
fclose(fid);
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function schnitzNum = findSchnitz(s,fr,cellnumber);
for schnitzNum = [1:length(s)]
  schnitz = s(schnitzNum);
  index = find( [schnitz.frames] == fr);
  if schnitz.cellno(index) == cellnumber
    return;
  end
end
schnitzNum = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
