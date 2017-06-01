function [problemCells, theOutputFilePath] = DJK_analyzeTracking(p, varargin);
% DJK_analyzeTracking analyses the complete tracking file. 
% **** if getting weird errors: 1) introduce if ~isempty tests (see also
%                                  NW2015-05 below)
%                               2) Try DJK_analyzeTracking_old
%
% Reports on :
% # cells tracked in complete lineage 
%   # cells in manualRange
%     * cells that are in segmentation but not in tree
%     * cells that are orphan (first frame always has orphans)
%     * cells that are barren (cell lineage ends before last frame)
%     * cells that move more than p.pixelsMoveDef pixels between frames 
%     * cells that grow more/less than p.pixelsLenDef pixels between frames 
%     * cells that change length more/less than 20 pixels after division
%                              and grow/shrink more than p.pixelAreaDiv
%
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
%                   Note: old length definition is used (fit of ellipse) ->
%                   majoraxislength. Using 'THIN' WOULD BE BETTER
%
%  'pixelsAreaDiv'   threshold change in px area at division. Extra criterion
%                  in addition to 20px length change because length change
%                  for bent cells (majoraxislength!) is a poor criterion
%                  and often false positive.
%                  default: 70
% p.noDisp          don't output what's written to screen..
%
%
% OUTPUT:
% problemcells; a matrix that contains: [schnitzes, frames, cellno] which appear
%               problematic. 
%               EDIT MW 2014/06/11: the problem occured for (frame) -> (frame+1)
%               Addition MW 2015/07: added cellno to this array
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

if ~isfield(p,'noDisp')
    p.noDisp = 0;
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
  for schnitzIdx = 1:length(schnitzcells)
    if (max(schnitzcells(schnitzIdx).frame_nrs) > fr_maximum) 
      fr_maximum = max(schnitzcells(schnitzIdx).frame_nrs);
    end
    if (min(schnitzcells(schnitzIdx).frame_nrs) < fr_minimum) 
      fr_minimum = min(schnitzcells(schnitzIdx).frame_nrs);
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

% If explicit pixelsLenDef is not given, use [-1 4] pixels
if ~existfield(p,'pixelsLenDef')
    p.pixelsLenDef = [-4 6];
end

% If explicit pixelsAreaDiv is not given, use 50 pixels
if ~existfield(p,'pixelsAreaDiv')
    p.pixelsAreaDiv = 70;
end

if ~isfield(p,'dontShowExtendedReport')
    p.dontShowExtendedReport = 0;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Open file to write results to
%--------------------------------------------------------------------------
theOutputFilePath = [p.DJK_saveDir p.movieName '-tracking.txt'];
fid = fopen(theOutputFilePath,'wt');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PREPARATION OF TRACKING DATA
%--------------------------------------------------------------------------
tic

dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
dispAndWrite(fid, [' * Saving results in : ' p.DJK_saveDir]);
dispAndWrite(fid, [' * Loaded tracking file from : ' p.lineageName]);

firstFrame = min(p.manualRange);
lastFrame = max(p.manualRange);

% Determine for each cell in schnitzcells whether really in tracking
% indicated by inTracking (schnitzcells(1).inTracking=true)
range = p.manualRange; % (MW 2014/06/11) N+1 fix
nrCellsInTracking = 0;
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).inTracking = false;
  if (intersect(range, schnitzcells(schnitzIdx).frame_nrs))
    schnitzcells(schnitzIdx).inTracking = true;
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
% MWc loop over schnitzes 
for schnitznum = 1:length(schnitzcells)
  % MWc loop over frames belonging to that schnitz
  for frame = schnitzcells(schnitznum).frame_nrs; % MW 2014/06/11 removal N+1 bug
    % MWc if frame within range
    if (frame<=lastFrame)
      % MWc increase count for that frame
      NrCellsInTrackingPerFrame(frame) = NrCellsInTrackingPerFrame(frame) + 1;
    end;
  end
end

disp(['Preparation done in ' num2str(toc) ' seconds']);

%% ------------------------------------------------------------------------
% ACTUAL CHECKS
%--------------------------------------------------------------------------

tic

% Find cells that are in segmentation but not in tree
segAndTrackDoNotMatch = zeros(lastFrame, 1);
for frameNum = p.manualRange
  if ~(NrCellsPerFrame(frameNum) == NrCellsInTrackingPerFrame(frameNum))
    segAndTrackDoNotMatch(frameNum) = 1;
  end
end

% Find orphans
nrOrphansInManualRange = 0;
problemSchnitzesOrphans = [];
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).orphan = true;
  if schnitzcells(schnitzIdx).P > 0;
    schnitzcells(schnitzIdx).orphan = false;
    problemSchnitzesOrphans(end+1) = schnitzIdx;
  else
    if schnitzcells(schnitzIdx).inTracking & schnitzcells(schnitzIdx).frame_nrs(1) > firstFrame % cells in first frame do not count % MW 2014/06/11 removal N+1 bug
      nrOrphansInManualRange = nrOrphansInManualRange + 1;
    end
  end
end
problemSchnitzesOrphans = unique(problemSchnitzesOrphans);

% output of problem cells (schnitz nr + frame)
problemCells = [];

% Find barren cells
framesWithBarrenCells = [];
problemSchnitzesBarren = [];
%nrBarrenInManualRange = 0;
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).barren = false;
  if schnitzcells(schnitzIdx).E==0 | schnitzcells(schnitzIdx).D==0;
    schnitzcells(schnitzIdx).barren = true;
    problemSchnitzesBarren(end+1) = schnitzIdx;
    if schnitzcells(schnitzIdx).inTracking & schnitzcells(schnitzIdx).frame_nrs(end) < lastFrame % cells in last frame do not count % MW 2014/06/11 removal N+1 bug
      framesWithBarrenCells = [framesWithBarrenCells schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
      problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
      % nrBarrenInManualRange = nrBarrenInManualRange + 1;
    end
  end
end
framesWithBarrenCells = unique(framesWithBarrenCells);
problemSchnitzesBarren = unique(problemSchnitzesBarren);

% Find cells moving > p.pixelsMoveDef pixels
framesWithCellsMoving = [];
problemSchnitzesMoving = [];
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).moving = [];
  cenx = schnitzcells(schnitzIdx).cenx_cent(1); % DJK 090410 cenx(i);
  ceny = schnitzcells(schnitzIdx).ceny_cent(1); % DJK 090410 ceny(i);
  for i = 2:length(schnitzcells(schnitzIdx).frame_nrs);
    if schnitzcells(schnitzIdx).inTracking % last frame should be taken in account (MW 2014/06/11)
      cenx_new = schnitzcells(schnitzIdx).cenx_cent(i); % DJK 090410 cenx(i);
      ceny_new = schnitzcells(schnitzIdx).ceny_cent(i); % DJK 090410 ceny(i);
      if sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 ) > p.pixelsMoveDef
        schnitzcells(schnitzIdx).moving = [schnitzcells(schnitzIdx).moving schnitzcells(schnitzIdx).frame_nrs(i-1)]; % MW 2014/06/11 here -1 needed
        problemSchnitzesMoving(end+1) = schnitzIdx;
        framesWithCellsMoving = [framesWithCellsMoving schnitzcells(schnitzIdx).frame_nrs(i-1)];        
        problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(i-1)];                   
%         disp([str3(cell) ' moved ' num2str(sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 )) ' pixels ']);
      end
      cenx = cenx_new; ceny = ceny_new;
    end
  end
end
framesWithCellsMoving = unique(framesWithCellsMoving);
problemSchnitzesMoving = unique(problemSchnitzesMoving);

% Find cells growing < > p.pixelsLenDef pixels
framesWithCellsGrowingTooLittle = [];
framesWithCellsGrowingTooMuch = [];
problemSchnitzesSlowGrowth = [];
problemSchnitzesFastGrowth = [];
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).growingTooLittle = [];
  schnitzcells(schnitzIdx).growingTooMuch = [];
  len = schnitzcells(schnitzIdx).len(1);
  for i = 2:length(schnitzcells(schnitzIdx).frame_nrs);
    if schnitzcells(schnitzIdx).inTracking % last frame should be taken in account (MW 2014/06/11)
      len_new = schnitzcells(schnitzIdx).len(i);
      if (len_new-len) < p.pixelsLenDef(1)
        schnitzcells(schnitzIdx).growingTooLittle = [schnitzcells(schnitzIdx).growingTooLittle schnitzcells(schnitzIdx).frame_nrs(i-1)]; % here -1 also required, MW 2014/06/11
        framesWithCellsGrowingTooLittle = [framesWithCellsGrowingTooLittle schnitzcells(schnitzIdx).frame_nrs(i-1)];
        problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(i-1)];
        problemSchnitzesSlowGrowth(end+1) = schnitzIdx;
      elseif (len_new-len) > p.pixelsLenDef(2) 
        schnitzcells(schnitzIdx).growingTooMuch = [schnitzcells(schnitzIdx).growingTooMuch schnitzcells(schnitzIdx).frame_nrs(i-1)];
        framesWithCellsGrowingTooMuch = [framesWithCellsGrowingTooMuch schnitzcells(schnitzIdx).frame_nrs(i-1)];
        problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(i-1)];
        problemSchnitzesFastGrowth(end+1) = schnitzIdx;
      end
      len = len_new;
    end
  end
end
framesWithCellsGrowingTooLittle = unique(framesWithCellsGrowingTooLittle);
framesWithCellsGrowingTooMuch = unique(framesWithCellsGrowingTooMuch);
framesWithCellsGrowingWeird = unique([framesWithCellsGrowingTooLittle framesWithCellsGrowingTooMuch]);
problemSchnitzesSlowGrowth = unique(problemSchnitzesSlowGrowth);
problemSchnitzesFastGrowth = unique(problemSchnitzesFastGrowth);

% Find cells that change length more/less than 20 pixels after division and
% have area change larger than p.pixelsAreaDiv
framesWithCellsChangingAfterDivision = [];
problemSchnitzesDivLenChange = [];
for schnitzIdx = 1:length(schnitzcells)
  schnitzcells(schnitzIdx).offspringToBig = false;
  D = schnitzcells(schnitzIdx).D;
  E = schnitzcells(schnitzIdx).E;
  if schnitzcells(schnitzIdx).inTracking & E>0 & D>0
    combined_length = schnitzcells(D).len(1) + schnitzcells(E).len(1);
    %length change
    if (schnitzcells(schnitzIdx).len(end) >  20 + schnitzcells(D).len(1) + schnitzcells(E).len(1)) | ...
       (schnitzcells(schnitzIdx).len(end) < -20 + schnitzcells(D).len(1) + schnitzcells(E).len(1))
        %area change (check only if schnitz contains this field -> not true for schnitz files older than 2013-12)
        if ~existfield(schnitzcells,'areapx')
            schnitzcells(schnitzIdx).offspringToBig = true;
            framesWithCellsChangingAfterDivision = [framesWithCellsChangingAfterDivision schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
            problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
            problemSchnitzesDivLenChange(end+1) = schnitzIdx;
        else
            AreaChange=abs(schnitzcells(schnitzIdx).areapx(end) -  ( schnitzcells(D).areapx(1) + schnitzcells(E).areapx(1)));
            if AreaChange>p.pixelsAreaDiv;
                schnitzcells(schnitzIdx).offspringToBig = true;
                framesWithCellsChangingAfterDivision = [framesWithCellsChangingAfterDivision schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
                problemCells = [problemCells ; schnitzIdx schnitzcells(schnitzIdx).frame_nrs(end)]; % MW 2014/06/11 removal N+1 bug
                problemSchnitzesDivLenChange(end+1) = schnitzIdx;
                %disp(['Parent schnitz: ' num2str(cell) ' Frame '   num2str(schnitzcells(cell).frame_nrs(end)-1)])
            end
        end
    end
  end
end
framesWithCellsChangingAfterDivision = unique(framesWithCellsChangingAfterDivision);
problemSchnitzesDivLenChange = unique(problemSchnitzesDivLenChange);

disp(['Done with analyzing problems in ' num2str(toc) ' seconds. Continuing output to file.']);

%--------------------------------------------------------------------------
% output of problem cells

problemCells = unique(problemCells, 'rows');

% MW addition, add a third column to the problemCells array, namely their
% cellnumbers, such that highlighting them becomes easier.
problemCells = padarray(problemCells,[0,1],'post');
for i = 1:numel(problemCells(:,1))
    schnitzIdx = problemCells(i,1);
    IdxInSchnitz = find(schnitzcells(schnitzIdx).frame_nrs == problemCells(i,2)); % label MW001
    problemCells(i,3) = schnitzcells(schnitzIdx).cellno(IdxInSchnitz);
end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% NOW DISPLAY ALL THESE PROBLEMS
%--------------------------------------------------------------------------
% get all framenumbers to write in which tracking step (e.g. 1->3 or
% 23->24) the problem occured
framesunique=unique([schnitzcells.frame_nrs]);
% The -1 below results in taking next frame and current frame as inspected
% frames, I think I solved this. (Variables are also more clearly named
% now.) - MW 2014/6/3
%framesunique=unique([schnitzcells.frame_nrs])-1;  %NW2013-12. Miraculeously, a shift of '-1' has to be introduced
                                               
% Display nr of cells in lineage file______________________________________
dispAndWrite(fid, ['-------------------------------------------------']);
if ~isempty(problemCells)
    dispAndWrite(fid, ['Frames with issues: ' mat2str(sort(unique(problemCells(:,2))))]);
    dispAndWrite(fid, ['Schnitzes with issues: ' mat2str(sort(unique(problemCells(:,1))))]);
end
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, 'Summary per category');
dispAndWrite(fid, ['problemSchnitzesBarren = [' num2str(problemSchnitzesBarren) ']']);
dispAndWrite(fid, ['problemSchnitzesMoving = [' num2str(problemSchnitzesMoving) ']']);
dispAndWrite(fid, ['problemSchnitzesSlowGrowth = [' num2str(problemSchnitzesSlowGrowth) ']']);
dispAndWrite(fid, ['problemSchnitzesFastGrowth = [' num2str(problemSchnitzesFastGrowth) ']']);
dispAndWrite(fid, ['problemSchnitzesDivLenChange = [' num2str(problemSchnitzesDivLenChange) ']']);
dispAndWrite(fid, '===');

dispAndWrite(fid, ['-------------------------------------------------']);

if p.dontShowExtendedReport
    
    dispAndWrite(fid, ['Not showing extended report. Activate by removing field p.dontShowExtendedReport.']);
    dispAndWrite(fid, ['Note that by pressing "s" in manual check seg your can also see schnitz colors & nrs..']);
    
else

    dispAndWrite(fid, ['Results:']);
    dispAndWrite(fid, ['--- ' str3(length(schnitzcells)) ' cells in tracking file']);
    dispAndWrite(fid, ['  | ']);
    dispAndWrite(fid, ['  |--- ' str3(nrCellsInTracking) ' of these are in manualRange']);
    dispAndWrite(fid, ['     | ']);

    % Display cells that are in segmentation but not in tree___________________
    dispAndWrite(fid, ['     |--- ' str3(sum(segAndTrackDoNotMatch)) ' frames lack cells in tracking, which are there in segmentation']);
    for frameNum = p.manualRange
      if segAndTrackDoNotMatch(frameNum)
        dispAndWrite(fid, ['     |  |--- frame ' str3(frameNum) ' has ' str3(NrCellsInTrackingPerFrame(frameNum)) ' / ' str3(NrCellsPerFrame(frameNum)) ' cells in tracking / segmentation']);
        for schnitzIdx = [1:NrCellsPerFrame(frameNum)]
          if ( findSchnitz(schnitzcells, frameNum, schnitzIdx) == 0) % 2014/06/11 MW +1 fix // == 0: edit other numbering issue          
              cellLabelNr = schnitzcells(schnitzIdx).cellno(find(schnitzcells(schnitzIdx).frame_nrs == frameNum)); % see label MW001 above
              dispAndWrite(fid, ['     |     |--- cellno ' str3(schnitzIdx) ' is missing in tracking (label=' num2str(cellLabelNr) ')' ]);
          end
        end
      end
    end

    % Display orphans (first frame always has orphans)_________________________
    dispAndWrite(fid, ['     |  ']);
    dispAndWrite(fid, ['     |--- ' str3(nrOrphansInManualRange) ' of cells are orphans']);
    for schnitzIdx = 1:length(schnitzcells)
      if schnitzcells(schnitzIdx).inTracking & schnitzcells(schnitzIdx).orphan & schnitzcells(schnitzIdx).frame_nrs(1) > firstFrame % cells in first frame do not count % MW 2014/06/11 removal N+1 bug
          cellLabelNr = schnitzcells(schnitzIdx).cellno(1); % see label MW001 above
            dispAndWrite(fid, ['     |  |--- schnitz ' str3(schnitzIdx) ' is an orphan and appears in frame ' str3(schnitzcells(schnitzIdx).frame_nrs(1)) '(label=' num2str(cellLabelNr) ')']); % MW 2014/06/11 removal N+1 bug
      end
    end

    % Display barren cells (cell lineage ends)_________________________________
    dispAndWrite(fid, ['     |  ']);
    % dispAndWrite(fid, ['     |--- ' str3(nrBarrenInManualRange) ' of barren cells in manualRange']);
    dispAndWrite(fid, ['     |--- ' str3(length(framesWithBarrenCells)) ' frames have barren cells']);
    if ~isempty(framesWithBarrenCells)
        for fr = framesWithBarrenCells
          out = ['     |  |--- in frame ' str3(fr) ' schnitz '];
          for schnitzIdx = 1:length(schnitzcells)
            if schnitzcells(schnitzIdx).inTracking & schnitzcells(schnitzIdx).barren & intersect(fr, schnitzcells(schnitzIdx).frame_nrs(end)) % MW 2014/06/11 removal N+1 bug

                % find label nr of cell (MW)
                cellLabelNr = schnitzcells(schnitzIdx).cellno(find(schnitzcells(schnitzIdx).frame_nrs == fr)); % see label MW001 above

              out = [out str3(schnitzIdx) ' (lbl=' num2str(cellLabelNr) ') '];
              %dispAndWrite(fid, ['     |  |--- schnitz ' str3(cell) ' is barren and disappears in frame ' str3(schnitzcells(cell).frame_nrs(end)-1)]);
            end
          end


          out = [out 'are barren'];
          dispAndWrite(fid, out);
        end
    end

    % Display moving cells (possibly wrong tracking)___________________________
    dispAndWrite(fid, ['     |  ']);
    dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsMoving)) ' frames have cells moving > ' num2str(p.pixelsMoveDef) ' pixels']);

    % Only loop when frames w. suspicious cells were detected
    if ~isempty(framesWithCellsMoving)

        % Loop over frames w. cells that seem to be moving (suspicious)    
        for fr = framesWithCellsMoving



          % output
          dispAndWrite(fid, ['     |  |--- in frame ' str3(fr) ' -> ' str3(fr+1)]); %NW2013-12 change display (fr)->(next fr)
          for schnitzIdx = 1:length(schnitzcells)

              % Get the index that corresponds to the schnitz with the frame of
              % interest, and then also find the next schnitz 
              % edit MW 2014/6/3: cleaner code. renaming.      
              %next_fr_idx =  find(schnitzcells(cell).frame_nrs==(fr+1)); % take +1 frame as reference
              %current_fr_idx = next_fr_idx - 1;      
              % TODO (2014/6/3): not entirely clear to me why we want next frame as a
              % reference, and not like this:
              current_fr_idx = find(schnitzcells(schnitzIdx).frame_nrs==(fr)); % MW 2014/06/11 let's do like this for now (TODO remove these comments.)
              next_fr_idx = current_fr_idx + 1;           

              if schnitzcells(schnitzIdx).inTracking & intersect(fr, schnitzcells(schnitzIdx).moving)                                   

                  % get x,y locations of centers of schnitzes
                  cenx = schnitzcells(schnitzIdx).cenx_cent(current_fr_idx); % MW 2014/6/3 renaming
                  ceny = schnitzcells(schnitzIdx).ceny_cent(current_fr_idx);
                  cenx_new = schnitzcells(schnitzIdx).cenx_cent(next_fr_idx);
                  ceny_new = schnitzcells(schnitzIdx).ceny_cent(next_fr_idx);

                  % find label nr of cell (MW)
                  cellLabelNr = schnitzcells(schnitzIdx).cellno(current_fr_idx); % see label MW001 above

                  % Pythagoras
                  distanceMoved = sqrt( (cenx_new-cenx)^2 + (ceny_new-ceny)^2 );
                  dispAndWrite(fid, ['     |  |  |--- schnitz ' str3(schnitzIdx) ' : ' num2str( round(distanceMoved) ) ' pixels (label=' num2str(cellLabelNr) ')']);

                  % MW 2014/6/3 || bugfix 2014/06/24
                  % missing frames might introduce extra movement, could tell user that
                  if ((find(schnitzcells(schnitzIdx).frame_nrs==next_fr_idx) - find(schnitzcells(schnitzIdx).frame_nrs==current_fr_idx)) ~= 1)
                      dispAndWrite(fid, ['     |  |--- (but missing frame detected, this might be the cause.)' ]);    
                  end

              end

          end

          dispAndWrite(fid, ['     |  |' ]);
        end
    end

    % Display weird growing cells (possibly wrong tracking / segmentation)_____
    dispAndWrite(fid, ['     |  ']);
    dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsGrowingWeird)) ' frames have cells growing < ' num2str(p.pixelsLenDef(1)) ' or > ' num2str(p.pixelsLenDef(2)) ' pixels']);
    if ~isempty(framesWithCellsGrowingWeird)
        for fr = framesWithCellsGrowingWeird       

          % output
          dispAndWrite(fid, ['     |  |--- in frame ' str3(fr) ' -> ' str3(fr+1)]);
          for schnitzIdx = 1:length(schnitzcells)

              % Get the index that corresponds to the schnitz with the frame of
              % interest, and then also find the next schnitz 
              % edit MW 2014/6/3: cleaner code. renaming.      
              %next_fr_idx =  find(schnitzcells(cell).frame_nrs==(fr+1)); % take +1 frame as reference
              %current_fr_idx = next_fr_idx - 1;      
              % TODO: not entirely clear to me why we want next frame as a
              % reference, and not like this:
              current_fr_idx = find(schnitzcells(schnitzIdx).frame_nrs==(fr)); % MW 2014/06/11
              next_fr_idx = current_fr_idx + 1;  % MW 2014/06/11

              % get indices schnitzes that grow too fast/too little
              idx = [find(schnitzcells(schnitzIdx).growingTooLittle==(fr))  find(schnitzcells(schnitzIdx).growingTooMuch==(fr))];

              if schnitzcells(schnitzIdx).inTracking & ~isempty(idx) % MW TODO 2014/6, probably this should be & ~isempty(idx)??

                % find label nr of cell (MW)
                cellLabelNr = schnitzcells(schnitzIdx).cellno(find(schnitzcells(schnitzIdx).frame_nrs == fr)); % see label MW001 above

                lengthIncrease = schnitzcells(schnitzIdx).len(next_fr_idx) - schnitzcells(schnitzIdx).len(current_fr_idx);
                dispAndWrite(fid, ['     |  |  |--- schnitz ' str3(schnitzIdx) ' : ' num2str( round(lengthIncrease) ) ' pixels (label=' num2str(cellLabelNr) ')']);

                % MW 2014/6/3 || bugfix 2014/06/24
                % missing frames might introduce extra movement, could tell user that
                if (find(schnitzcells(schnitzIdx).frame_nrs==next_fr_idx)-find(schnitzcells(schnitzIdx).frame_nrs==current_fr_idx) ~= 1)
                    dispAndWrite(fid, ['     |  |--- (but missing frame detected, this might be the cause.)' ]);    
                end          

              end

          end

          dispAndWrite(fid, ['     |  |' ]);
        end
    end

    % Display moving whose length change after division (possibly wrong tracking)
    dispAndWrite(fid, ['     |  ']);
    dispAndWrite(fid, ['     |--- ' str3(length(framesWithCellsChangingAfterDivision)) ' frames have cells with length change <> 20 pixels after division']);

    % loop over frames when there are frames that contain suspicious cells
    if ~isempty(framesWithCellsChangingAfterDivision)
        for fr = framesWithCellsChangingAfterDivision

          % get next frame nr and current frame nr (MW 2015/07 fix)
          current_fr_idx = fr;
          tableIdx = find(framesunique==fr)+1;
          if tableIdx>numel(framesunique)
                warning(['The strangest thing happened, skipping remark about ' num2str(current_fr_idx) ' since it''s last frame']); 
                continue
          end;
          next_fr_idx = framesunique(tableIdx);

          % output
          dispAndWrite(fid, ['     |  |--- in frame ' str3(current_fr_idx) ' -> ' str3(next_fr_idx)]); 
          for schnitzIdx = 1:length(schnitzcells)

            % determine whether this schnitz is suspicious
            if schnitzcells(schnitzIdx).inTracking & schnitzcells(schnitzIdx).offspringToBig & fr==schnitzcells(schnitzIdx).frame_nrs(end) % MW 2014/06/11 removal N+1 bug

              % Obtain children
              D = schnitzcells(schnitzIdx).D;
              E = schnitzcells(schnitzIdx).E;

              % Calculate length increase
              lengthIncrease = schnitzcells(D).len(1) +  schnitzcells(E).len(1) - schnitzcells(schnitzIdx).len(end);          

              % find label nr of cell (MW)
              cellLabelNr = schnitzcells(schnitzIdx).cellno(find(schnitzcells(schnitzIdx).frame_nrs == current_fr_idx)); % see label MW001 above

              % Output
              dispAndWrite(fid, ['        |  |--- parent schnitz ' str3(schnitzIdx) ' : ' num2str( round(lengthIncrease) ) ' pixels (label=' num2str(cellLabelNr) ')']);

              %{
              % MW 2014/6/3 || bugfix 2014/6/24
              % missing frames might introduce extra movement, could tell user that
              if (find(schnitzcells(cell).frame_nrs==next_fr_idx)-find(schnitzcells(cell).frame_nrs==current_fr_idx) ~= 1)
                  dispAndWrite(fid, ['     |  |--- (but missing frame detected, this might be the cause.)' ]);    
              end          
              %}

            end        

          end      
          dispAndWrite(fid, ['        |' ]);
        end
    end
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
% Find schnitznumber associated with certain framenumber and cellnumber.
% Find schnitznumber that contains a certain cell (identified by its
% number cellno), in a certain frame (identified by its number fr).

% loop over schnitzes
for schnitzNum = [1:length(s)]
  % get schnitz
  schnitz = s(schnitzNum);
  
  % get idx corresponding to fr for schnitzarray
  index = find( [schnitz.frame_nrs] == fr);
  
  % check whether current schnitz is the cell we're looking for
  if schnitz.cellno(index) == cellnumber
    return;
  end
end

% schnitzNum = -1; % this makes no sense - MW 2014/06/11 (edit also with
% usage).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








