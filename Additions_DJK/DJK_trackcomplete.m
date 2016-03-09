function p = DJK_trackcomplete (p, varargin)
% DJK_TRACKCOMPLETE     Track cells across frames in one complete movie.
% 
%   Same as TRACKCOMPLETE but allows SingleCell track method. To use, give
%   argument 'trackMethod','singleCell'
%  
%   TRACKCOMPLETE allows users to track cells identified by cell segmentation 
%   across all frames for a movie.  Cells in one frame are matched to cells 
%   in adjacent frames, and cell division is recognized by associating each 
%   parent cell with its children.
%   
%   TRACKCOMPLETE(P,'Field1',Value1,'Field2',Value2,...) also performs cell
%   tracking using segmentation results, but the extra arguments permit users 
%   to adjust any parameters describing the movie or parameters controlling 
%   the cell tracking process.  The extra arguments can overwrite any specific 
%   parameters provided in P by setting P.Field1 = Value1, P.Field2 = Value2, 
%   etc.  Thus any/all schnitzcells parameter values can be defined in the 
%   function call via these optional field/value pairs.  (This is in the style 
%   of setting MATLAB properties using optional property/value pairs.)  
%   
%   TRACKCOMPLETE produces an output file (a MATLAB binary file) that 
%   contains the 'schnitzcells' variable describing each cell's lineage. The 
%   path name of this file defaults to [p.tracksDir,p.movieName,'_lin.mat'], 
%   but can be specified by providing an optional 'lineageName' field in the 
%   schnitcells parameter structure, or by providing an extra 
%   'lineageName',value pair of arguments to this function. 
%   
%   TRACKCOMPLETE returns a struct (1x1 struct array) referred to as the 
%   schnitzcells parameter structure that contains fields and values 
%   contained in P, including unchanged/original parameters plus any of those 
%   added or overridden in the list of properties & values provided in the 
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control trackcomplete
% 
%   lineageName   the schnitzcells cell lineage structure is written to this 
%                 file, by default named [p.trackDir p.moveiName '_lin.mat']
%   
%   trackMethod   name of tracking method, one of {'tps','original'}; tps is
%                 the default tracking method
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
%   transMax      maximum translation accepted; default is 30
%   
%   len           size of subimage used for cross correlation; default is 80
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
    fieldName = DJK_schnitzfield(varargin{i}); % DJK 090210
    p.(fieldName) = varargin{i+1};
  end
end

% handle tracking method- switch control over to new TPS tracker if needed
if ~existfield(p,'trackMethod')
  % TPS is now the default tracker as of 2005-06-26
  p.trackMethod = 'tps';
end

if strcmp(upper(p.trackMethod),'TPS') 
  % make sure tps tracker is found in user's path
  if isempty(which('tracker_tps'))
    error(['TPS tracker path not found! '...
           'Add schnitzcells tracking_tps folder to your path.'])
  end
  % if using TPS tracking, switch control over to that module
  p = DJK_tracker_tps(p); % DJK 071129
  return;
end

if strcmp(upper(p.trackMethod),'SINGLECELL') 
  % if using single Cell tracking, switch control over to that module
  p = DJK_tracker_singleCell(p);
  return;
end

if ~strcmp(upper(p.trackMethod),'ORIGINAL')
  error(['Unknown trackMethod ' p.trackMethod ' ; ' ...
         'Please use ''original'' or ''tps'''])
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

if ~existfield(p,'transMax')
  p.transMax = 30;
end
if ~existfield(p,'len')
  p.len = 80;
end

%% Function appears to be unused - MW 2014/06/24
%{

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

% Internally, code is set up to use the track range in reverse
numFrames = max(size(p.trackRange));

% Short cut: If segmentation is fillinator output, we can take a shortcut
% Note: Now that 'original' is not the default, fillinator users have to 
% explicitly specify trackMethod original
% Note: We should make this it's own alternate trackMethod
if existfield(p,'fillinatorFile')
    disp('Checking if this segmentation is fillinator output...')

    % create schnitzcells fields assuming it *is* a fillinator segmentation
    fillinatorSegmentation = true;
    schnitzcells.P = -1;
    schnitzcells.frames = [];
    schnitzcells.cellno = [];
    schnitzcells.cenx = [];
    schnitzcells.ceny = [];
    schnitzcells.D = -1;
    schnitzcells.E = -1;
    schnitzcells.N = numFrames;
    for frameNum = p.trackRange
        todayname= [p.segmentationDir, outprefix, str3(frameNum), '.mat'];
        % modified nitzan 2005June24
        % clear Lc;
        % load(todayname);
        clear Lc LNsub;
        load(todayname);
        if exist('Lc')~=1 & p.trackUnCheckedFrames   
            Lc = LNsub;
        end
        % end nitzan 2005June24
        MaxCell = max2(Lc);
        if MaxCell == 0
            error(['Segmentation for frame ',num2str(frameNum),...
                   'has no cells!  Go fix this before tracking!']);
        end
        if MaxCell > 1
            disp('The segmentation is NOT fillinator output, proceeding.')
            fillinatorSegmentation = false;
            break;
        end
        % up till now looks like fillinator output, so keep gathering info
        r = regionprops(Lc,'Centroid');
        schnitzcells.cenx = [schnitzcells.cenx r(1).Centroid(1)];
        schnitzcells.ceny = [schnitzcells.ceny r(1).Centroid(2)];
        schnitzcells.cellno = [schnitzcells.cellno 1];
        schnitzcells.frames = [schnitzcells.frames frameNum];
    end
    if fillinatorSegmentation
        disp('Writing schnitzcell lineage for fillinator segmentation...')
        schnitzcells.frames = schnitzcells.frames + 1; % JCR: Hack- schnitzedit
        save(p.lineageName,'schnitzcells');
        return;
    end
end

partialname= [p.tracksDir, 'partial\', 'partialtracks.mat'];
usepartial = false;


% NOTE: we don't include the first trackRange frame number, will not enter 
% tracking loop for this first frame number since we only process frame *pairs*
reverseTrackRange = p.trackRange(numFrames:-1:2);  

% % Load partial tracking results if they exist and if the trackrange matches
% partialname= [p.partialDir, 'partialtracks.mat'];
% usepartial = false;
% if exist(partialname)==2 & ~p.overwrite % instead of: "if exist(partialname)==2"
%     partialtracks = load(partialname);
%     disp('Checking partial tracks file...');
%     if (length(partialtracks.trackRange) == length(p.trackRange))
%         if (partialtracks.trackRange == p.trackRange)
%             load(partialname);
%             usepartial = true;
%             disp(['Partial tracks file matches, restarting tracking at ' ...
%                   num2str(lastFrameNumTracked)]);
%         end
%     else
%         disp('Partial tracks file track range does not match; ignoring');
%     end
% end

% CARRY OUT TRACKING
i = length(p.trackRange);  % i = row index into cellmat, NOT frameNum
for frameNum = reverseTrackRange
    todayname= [p.segmentationDir, outprefix, str3(frameNum), '.mat'];
    mynum_t= str3(frameNum);
    yesterdayFrameNum = p.trackRange(find(p.trackRange==frameNum)-1);
    yesterdayname = [p.segmentationDir, outprefix, str3(yesterdayFrameNum)];
    mynum_y= str3(yesterdayFrameNum);

    if (usepartial)
        if (frameNum >= lastFrameNumTracked)
            disp(['Skipping previously tracked frame pair ',...
		    mynum_y,' , ',mynum_t,'  (use p.overwrite=1 to redo)']);
	    i = i-1;
	    continue
        end
    end
    disp(['Starting to process frame pair ',mynum_y,' , ',mynum_t]);

    % current image
    clear LNsub Lc;
    load(todayname);
    % modified nitzan 2005June24
    % today= zeros([1024 1280]); % replaced by the following:
    if exist('phaseFullSize')==1
        today= zeros(phaseFullSize);
    elseif isfield(p,'fullsize');
        today= zeros(p.fullsize);
    else
        today= zeros([1024 1280]);
        disp('Using default image size.');
    end
    % end nitzan 2005June24
    
    % nitzan 2005June24 - i moved this line into the if scope below
    % rect(1)=1;rect(2)=1;rect(3)=size(Lc,1);rect(4)=size(Lc,2);
    
    if exist('Lc')
        rect(1)=1;rect(2)=1;rect(3)=size(Lc,1);rect(4)=size(Lc,2); % nitzan2005june24
        today(rect(1):rect(3), rect(2):rect(4))= imerode(Lc, se);
    else      
        rect(1)=1;rect(2)=1;rect(3)=size(LNsub,1);rect(4)=size(LNsub,2); % nitzan2005june24
        today(rect(1):rect(3), rect(2):rect(4))= imerode(LNsub, se);
    end;

    % previous image
    clear LNsub Lc;
    load(yesterdayname);
    % modified nitzan 2005June24
    % yesterday= zeros([1024 1280]); % replaced by the following:
    if exist('phaseFullSize')==1
        yesterday= zeros(phaseFullSize);
    elseif isfield(p,'fullsize');
        yesterday= zeros(p.fullsize);
    else
        yesterday= zeros([1024 1280]);
        disp('Using default image size.');
    end
    % end nitzan 2005June24
    
    % nitzan 2005June24 - i moved this line into the if scope below
    % rect(1)=1;rect(2)=1;rect(3)=size(Lc,1);rect(4)=size(Lc,2);
    
    if exist('Lc')
        rect(1)=1;rect(2)=1;rect(3)=size(Lc,1);rect(4)=size(Lc,2); % nitzan2005june24
        yesterday(rect(1):rect(3), rect(2):rect(4))= imerode(Lc, se);
    else
        rect(1)=1;rect(2)=1;rect(3)=size(LNsub,1);rect(4)=size(LNsub,2); % nitzan2005june24
        yesterday(rect(1):rect(3), rect(2):rect(4))= imerode(LNsub, se);
    end;

    % find centre of mass transformation
    [cmtX, cmtY]= centremass(today);
    [cmyX, cmyY]= centremass(yesterday);    
    trans= round([cmyX - cmtX, cmyY - cmtY]);

    % extend images so that trans doesn't move coordinates off image
    padsize= round(0.5*(max(abs(trans)) + p.len));
    % today
    todaypad= zeros(size(today) + 2*padsize);
    todaypad(padsize + 1:end - padsize, padsize + 1: end - padsize)= today;
    % yesterday
    yesterdaypad= zeros(size(today) + 2*padsize);
    yesterdaypad(padsize + 1:end - padsize, padsize + 1: end - padsize)= yesterday;
    ycells= unique(yesterday(:));
    ycells(1)= [];

    % TRACK EACH CELL
    % backtracking from today to yesterday
    disp(['backtracking cells: ',num2str(length(tcells))]);   
    clear bt scb;
    for j= 1:length(tcells)
%         disp(['  bcell no ',num2str(j)]);
        if tcells(j) > 0           
            % store prediction in cellmat
            [bt(j), scb(j)]= backtrack(today, yesterdaypad, tcells(j), trans, p.transMax, p.len, padsize);
        else
            % a zero is tracked to a zero
            bt(j)= 0;
            scb(j)= -1;
        end;
    end;

    % forwardtracking from yesterday to today
    disp(['forwardtracking cells: ',num2str(length(ycells))]);   
    clear ft scf;
    for j= 1:length(ycells)
%         disp(['  fcell no ',num2str(j)]);
        if ycells(j) > 0           
            % store prediction in cellmat
            [ft{j}, scf(j)]= forwardtrack(todaypad, yesterday, ycells(j), trans, p.transMax, p.len/2, padsize);
        else
            % a zero is tracked to a zero
            ft{j}= 0;
            scf(j)= -1;
        end;
    end;

    % CHECK TRACKING
    btold= bt;
    % find orphan cells
    %     orphans= setdiff(1:max2(yesterday), unique(bt));
    uycells = setdiff(unique(yesterday),0);
    orphans = setdiff(uycells, unique(bt));

    orphansold{i}= orphans;
    % reinstate orphans using forwardtracking
    odel= [];
    toomany{i}= [];
    noorphans= length(orphans);
    currorphan= 1;
    while noorphans > 0
        disp([' ',num2str(currorphan),': checking orphan ',num2str(orphans(currorphan))]);
        loc= find(orphans(currorphan) == ycells);
        parentcell= ft{loc};
        % include all neighbouring cells of forwardtrack predictions
        if length(parentcell) == 1
            parentcell= [parentcell neighbouringcells(today, parentcell)];
        else
            parentcell= [parentcell unique([neighbouringcells(today, parentcell(1)) ...
                        neighbouringcells(today, parentcell(2))])];
        end;

        % run through potential parents checking if any can be backtracked to orphan
        for k= 1:length(parentcell)
            btt= [];
            if (parentcell(k)==0)
                disp(' ERROR: TRACKCOMPLETE CANT BACKTRACK BACKGROUND! (continuing)');
                continue
            end
            btt(1)= backtrack(today, yesterdaypad, parentcell(k), trans, p.transMax, p.len/4, padsize);
            btt(2)= backtrack(today, yesterdaypad, parentcell(k), trans, p.transMax, p.len/2, padsize);
            btt(3)= backtrack(today, yesterdaypad, parentcell(k), trans, p.transMax, 2*p.len, padsize);
            if ismember(orphans(currorphan), unique(btt)) & parentcell(k)>0
                oldprediction= bt(parentcell(k));
                % replace old prediction with orphan
                bt(parentcell(k))= orphans(currorphan);
                ft{loc}= parentcell(k);
                odel= [odel currorphan];
                % check if replacement of old prediction results in new orphan 
                noappearances= find(bt == oldprediction); 
                if length(noappearances) == 1    
                    if ~ismember(oldprediction, orphans)  
                        disp([' adding new orphan ', num2str(oldprediction)]);
                        orphans= [orphans oldprediction];
                        noorphans= noorphans + 1;
                    end;
                end;
                % jump out of parentcell loop once orphan has been parented
                break;
            end;
        end;     

        currorphan= currorphan + 1;
        noorphans= noorphans - 1;
    end;

    % delete orphans successfully parented
    odel= unique(odel);
    orphans(odel)= [];
    % store any remaining orphans
    if ~isempty(orphans)
        disp([' found ',num2str(length(orphans)),' orphans ',num2str(orphans)]);
        disp([' frame= ',num2str(frameNum)]);
        orphansnew{i}= orphans;
    else
        orphansnew{i}= [];
    end;

    % find cells that aren't reciprocally tracked
    badcells{i}= []; 
    for j= 1:length(tcells)
        if ~ismember(tcells(j), ft{find(bt(j) == ycells)})
            badcells{i}= [badcells{i} tcells(j)];
        end;
    end;  

    % replace badcell with neighbouring cell which is correctly forward tracked
    % only replace badcell if centroid of backcell was initially tracked to empty space
    bdel= [];
    for j= 1:length(badcells{i})
        disp([' checking badcell ', num2str(j)]);
        locb= badcells{i}(j);
        if scb(locb) == 1 
            ns= neighbouringcells(yesterday, bt(locb));
            for k= 1:length(ns)
                locf= find(ns(k) == ycells);
                if ismember(badcells{i}(j), ft{locf}) & scf(locf) == 0
                    bt(locb)= ns(k);
                    bdel= [bdel j];
                end;
            end;
        end;
    end;
    badcells{i}(bdel)= [];

    if i > 1
        % INSERT INTO CELLMAT      
        for j= 1:length(tcells)
            cellmat(i-1, find(cellmat(i, :) == j))= bt(j);
        end;
        % RESET TCELLS
        tcells= unique(cellmat(i-1, :));
        % nitzan's addition July 18th - lines 214-215.
        tcells(tcells==0)=[];
    end;

    % CHECK FOR ANY REMAINING ORPHANS OR ANY NEW ORPHANS CREATED
    clear orphans;
    orphans= setdiff(1:max2(yesterday), tcells);
    if ~isempty(orphans)
        disp(['final number of orphans ',num2str(length(orphans))]);    
        % store orphans
        n1= size(cellmat,2);
        cellmat(i:end, n1+1:n1+length(orphans))= 0;
        cellmat(i-1, n1+1:n1+length(orphans))= orphans;
        % find new tcells
        tcells= unique(cellmat(i-1, :));
        % nitzan's addition July 18th - lines 214-215.
        tcells(tcells==0)=[];
    end;

    lastFrameNumTracked = frameNum;
    trackRange = p.trackRange;
    save(partialname, 'lastFrameNumTracked', 'trackRange', 'cellmat', 'badcells', 'orphansold', 'tcells');
    
    i = i-1;  % index into cellmat, NOT frame number
end;

% RECORD DATA
trackRange = p.trackRange;
save(p.trackName, 'trackRange', 'cellmat', 'badcells', 'orphansold');

% Convert the tracking results into a schnitzcells structure
p = makeschnitz(p);
%}
