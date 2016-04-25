function p = DJK_tracker_djk (p, varargin)
% DJK_tracker_djk  Track cells across frames in one complete movie in original fashion.
% 
%   Copied from DJK_TRACKCOMPLETE:
%     * taken tps out
%     * 
%
%--------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control trackcomplete
% 
%   lineageName   the schnitzcells cell lineage structure is written to this 
%                 file, by default named [p.trackDir p.moveiName '_lin.mat']
%   
%   manualRange   range of frame numbers to track across; by default all 
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
%
%  [NW 2013-09: len and transMax are non-existent parameters in the code!
%  There is also nothing like a cross correlation but the tracking
%  connectionas are determined by least square distance between the
%  thin-line of new and old cells]
%
%   transMax      maximum translation accepted; default is 30
%   
%   len           size of subimage used for cross correlation; default is 80
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1;
functionName = 'DJK_tracker_djk';

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
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
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
  disp('Please use p.overwrite instead of p.override.');
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
  disp('Set manual range to all available data.');
end

% Keep only the frames in the range that contain a corrected segmentation (unless we're tracking un-checked frames)
manualRangeChecked = [];
for frameNum = p.manualRange
  clear Lc;
  load([p.segmentationDir,p.movieName,'seg',str3(frameNum)]);
  if exist('Lc')==1 | p.trackUnCheckedFrames 
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


%--------------------------------------------------------------------------
% CARRY OUT TRACKING -> LOOP OVER FRAME PAIRS
%--------------------------------------------------------------------------
frameNum = 0;
dataLeft_previousRound = false;
for count = 2:length(p.manualRange);
  frameNum_previousRound = frameNum;
  
  % today is frameNum
  frameNum = p.manualRange(count);

  % yesterday is frameNum-1, will track between these frames
  yesterdayFrameNum = p.manualRange(count-1);

  % output for each frame to frame track, will be saved here
  trackOutputFile = [p.tracksDir,p.movieName,'-djk-output-',str3(yesterdayFrameNum),'-to-',str3(frameNum),'.txt'];

  % segmentation files used for tracking
  todaySegFile = [p.segmentationDir, outprefix, str3(frameNum), '.mat'];
  yesterdaySegFile = [p.segmentationDir, outprefix, str3(yesterdayFrameNum), '.mat'];
  
  
  %------------------------------------------------------------------------
  % CHECK WHETHER PAIR NEEDS TO BE TRACKED
  %------------------------------------------------------------------------
  % assume this pair of frames still need to be tracked
  needToTrack = true;

  if exist(trackOutputFile)==2 & ~p.overwrite
    % if trackOutputFile already exists, might not need to track
    needToTrack = false;
    
    % if trackOutputFile is younger than segmentation files, means that segmentation has been updated, and we want to track again
    info_yesterdaySegFile = dir(yesterdaySegFile);
    info_todaySegFile = dir(todaySegFile);
    info_trackOutputFile = dir(trackOutputFile);
    datenumber_yesterdaySegFile = datenum(info_yesterdaySegFile.date);
    datenumber_todaySegFile = datenum(info_todaySegFile.date);
    datenumber_trackOutputFile = datenum(info_trackOutputFile.date);
    if datenumber_yesterdaySegFile>datenumber_trackOutputFile | datenumber_todaySegFile>datenumber_trackOutputFile
      needToTrack = true;
    end
  end
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % GET DATA FOR TRACKING
  %------------------------------------------------------------------------
  % load centered and original full size version of each Lc
  [Lc_today_fullsize_centered, Lc_today_fullsize, Lc_today]             = loadLcData(todaySegFile); % loadCenteredSegData2(todaySegFile);
  [Lc_yesterday_fullsize_centered, Lc_yesterday_fullsize, Lc_yesterday] = loadLcData(yesterdaySegFile); % loadCenteredSegData2(yesterdaySegFile);
  
  % TODO: it would better to save cenx & ceny, from thin
  
  % Need to save some properties of cells for recalc_schnitz lateron
  rp  = regionprops(Lc_today,'Centroid','Orientation','MajorAxisLength');
  rp2 = regionprops(Lc_today_fullsize_centered,'Centroid');
  num_pts = size(rp,1);
  for j=1:num_pts
    Points(j).cenx      = rp(j).Centroid(1);
    Points(j).ceny      = rp(j).Centroid(2);
    Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
    Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
    Points(j).ang       = rp(j).Orientation;
    Points(j).len       = rp(j).MajorAxisLength;
    Points(j).cellno    = j;
  end
  opts{count}=Points(1:num_pts);
  
  % in case of first run, one extra
  if count == 2
    rp  = regionprops(Lc_yesterday,'Centroid','Orientation','MajorAxisLength');
    rp2 = regionprops(Lc_yesterday_fullsize_centered,'Centroid');
    num_pts = size(rp,1);
    for j=1:num_pts
      Points(j).cenx    = rp(j).Centroid(1);
      Points(j).ceny    = rp(j).Centroid(2);
      Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
      Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
      Points(j).ang     = rp(j).Orientation;
      Points(j).len     = rp(j).MajorAxisLength;
      Points(j).cellno  = j;
    end
    opts{count-1}=Points(1:num_pts);
  end
  
  % tell what is happening
  fprintf([' * frame pair ',str3(yesterdayFrameNum),' -> ', str3(frameNum) ' : ']);
  if ~needToTrack
    dataLeft_previousRound = false;
    fprintf(1,' -> Skipping, cause seg older than previous tracking (use p.overwrite=1 to redo)\n');
    continue
  end
  
  % cell numbers in each images (background is also in here with value 0)
  cellnos_today      = unique(Lc_today_fullsize_centered);
  cellnos_yesterday  = unique(Lc_yesterday_fullsize_centered);
  fprintf(1,'.');

  % determine 3 coordinates for each cellno: 1/4, 1/2 & 3/4 of thin
  % get coordinates of yesterday
  coordinates_yesterday = [];
  if yesterdayFrameNum==frameNum_previousRound & dataLeft_previousRound
    % can use old data
    coordinates_yesterday = coordinates_today;
  else
    for cellno = cellnos_yesterday'
      if cellno==0, continue; end % background
      coordinates_yesterday(cellno,:) = getCentroidsOfCell(Lc_yesterday_fullsize_centered, cellno);
    end
    coordinates_yesterday = normalizeCentroids(coordinates_yesterday, Lc_yesterday_fullsize_centered);
  end
  fprintf(1,'.');

  % get coordinates of today
  coordinates_today = [];
  for cellno = cellnos_today'
    if cellno==0, continue; end % background
    coordinates_today(cellno,:) = getCentroidsOfCell(Lc_today_fullsize_centered, cellno);
  end
  coordinates_today = normalizeCentroids(coordinates_today, Lc_today_fullsize_centered);
  dataLeft_previousRound = true;
  fprintf(1,'.');

  % DEBUGGING: show image
  if false
    image_today = zeros(size(Lc_today_fullsize_centered));
    for i = 1:size(coordinates_today,1)
      centroids = coordinates_today(i,:);
      image_today( round(centroids(2)) , round(centroids(3)) ) = centroids(1);
      image_today( round(centroids(4)) , round(centroids(5)) ) = centroids(1);
      image_today( round(centroids(6)) , round(centroids(7)) ) = centroids(1);
    end
    figure(11); PN_imshowlabel(p,image_today,[],[],[]);

    image_yesterday = zeros(size(Lc_yesterday_fullsize_centered));
    for i = 1:size(coordinates_yesterday,1)
      centroids = coordinates_yesterday(i,:);
      image_yesterday( round(centroids(2)) , round(centroids(3)) ) = centroids(1);
      image_yesterday( round(centroids(4)) , round(centroids(5)) ) = centroids(1);
      image_yesterday( round(centroids(6)) , round(centroids(7)) ) = centroids(1);
    end
    figure(12); PN_imshowlabel(p,image_yesterday,[],[],[]);
    
%     figure(11); DJK_imshowlabel(Lc_today_fullsize);
%     figure(12); DJK_imshowlabel(Lc_yesterday_fullsize);
    figure(13); PN_imshowlabel(p,Lc_today_fullsize_centered,[],[],[]);
    figure(14); PN_imshowlabel(p,Lc_yesterday_fullsize_centered,[],[],[]);
    
    pause; close(11); close(12); close(13); close(14);
    DJK_writeSegImage(image_today,'image_today.png');
    DJK_writeSegImage(image_yesterday,'image_yesterday.png');
%     DJK_writeSegImage(Lc_today_fullsize,'Lc_today_fullsize.png');
%     DJK_writeSegImage(Lc_yesterday_fullsize,'Lc_yesterday_fullsize.png');
    DJK_writeSegImage(Lc_today_fullsize_centered,'Lc_today_fullsize_centered.png');
    DJK_writeSegImage(Lc_yesterday_fullsize_centered,'Lc_yesterday_fullsize_centered.png');
  end
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % CARRY OUT TRACKING
  %------------------------------------------------------------------------
  trackLink = [];
  
  % Now for each today coordinate, find closest yesterday coordinate
  yesterday_cellno = [coordinates_yesterday(:,1); coordinates_yesterday(:,1); coordinates_yesterday(:,1)];
  yesterday_x = [coordinates_yesterday(:,2); coordinates_yesterday(:,4); coordinates_yesterday(:,6)]; 
  yesterday_y = [coordinates_yesterday(:,3); coordinates_yesterday(:,5); coordinates_yesterday(:,7)]; 
  for cellno = cellnos_today'
    if cellno==0, continue; end % background
    cellno_x = coordinates_today(cellno,2); cellno_y = coordinates_today(cellno,3);
    dist_square = (yesterday_x-cellno_x).*(yesterday_x-cellno_x) + (yesterday_y-cellno_y).*(yesterday_y-cellno_y);
    [trash,index] = sort(dist_square);
    trackLink = [trackLink; yesterday_cellno(index(1)) 0 0 cellno];
  end
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % CORRECT ERRORS IN TRACKING: BARREN CELLS
  %------------------------------------------------------------------------
  % now loop over all cellno in yesterday, if any does not occur, correct
  for cellno = cellnos_yesterday'
    if cellno==0, continue; end % background
    
    % how often does this cellno come back?
    cellno_i = find(trackLink(:,1)==cellno);
    nr_occur = length(cellno_i);
    
    % if cellno occurs 0 times, link to closest dividing cell
    if nr_occur==0
      % loop over cells in yesterday that are linked to >= 2 cells
      today_dividing_list = [];
      for c = cellnos_yesterday'
        % get the cells they link to
        linked_cells = find(trackLink(:,1)==c);
        if length(linked_cells) >= 2
          % now add the cells they link to today_dividing_list
          today_dividing_list = [today_dividing_list linked_cells'];
        end
      end
      
      if length(today_dividing_list) > 0
        % calc distances and pick closest
        cellno_x = coordinates_yesterday(cellno,2); cellno_y = coordinates_yesterday(cellno,3);
        today_dividing_cellno = [coordinates_today(today_dividing_list,1); coordinates_today(today_dividing_list,1); coordinates_today(today_dividing_list,1)];
        today_dividing_x      = [coordinates_today(today_dividing_list,2); coordinates_today(today_dividing_list,4); coordinates_today(today_dividing_list,6)]; 
        today_dividing_y      = [coordinates_today(today_dividing_list,3); coordinates_today(today_dividing_list,5); coordinates_today(today_dividing_list,7)]; 
        dist_square = (today_dividing_x-cellno_x).*(today_dividing_x-cellno_x) + (today_dividing_y-cellno_y).*(today_dividing_y-cellno_y);
        [trash,index] = sort(dist_square);

        fprintf(['\n                               ' str3(cellno) ' in fr' str3(yesterdayFrameNum) ' was not linked. Corrected: ' num2str(trackLink(today_dividing_cellno(index(1)),:)) ' to ' num2str([cellno 0 0 today_dividing_cellno(index(1))])]);
        % update closest
        trackLink(today_dividing_cellno(index(1)),:) = [cellno 0 0 today_dividing_cellno(index(1))];
      else
        fprintf(['\n                               ' str3(cellno) ' in fr' str3(yesterdayFrameNum) ' was not linked. Could not be corrected, cause no dividing cells!']);
      end
    end
  end
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % CORRECT ERRORS IN TRACKING: MULTIPLE LINKS
  %------------------------------------------------------------------------
  % now loop over all cellno in yesterday, if any > twice, correct
  for cellno = cellnos_yesterday'
    if cellno==0, continue; end % background
    
    % how often does this cellno come back?
    cellno_i = find(trackLink(:,1)==cellno);
    nr_occur = length(cellno_i);
    
    % if cellno occurs more than 2 times, it must be linked to extra cell(s) -> correct for this
    if nr_occur>2 
      fprintf(['\n                               ' str3(cellno) ' in fr' str3(yesterdayFrameNum) ' is linked to ' num2str(nr_occur) ' cells in the next frame (' str3(frameNum) ')']);
    end
    
    % do correction, untill nr_occur==2
    while nr_occur>2 
      % get list of cells in yesterday frame which are linked to no cells
      barren_list = [];
      for cellno_yes = cellnos_yesterday'
        if cellno_yes==0, continue; end % background
        if length(find(trackLink(:,1)==cellno_yes)) == 0
          barren_list = [barren_list cellno_yes];
        end
      end
      % if no more barren cells, will continue with cells linked to 1 other cell
      if length(barren_list) == 0
        for cellno_yes = cellnos_yesterday'
          if cellno_yes==0, continue; end % background
          if length(find(trackLink(:,1)==cellno_yes)) == 1
            barren_list = [barren_list cellno_yes];
          end
        end
      end
        
      % calc distances and pick closest
      closest_distance          = realmax('double');
      closest_cellno_yesterday  = 1;
      closest_cellno_today      = 1;
      for i = 1:nr_occur
        cellno_x = coordinates_today(cellno_i(i),2); cellno_y = coordinates_today(cellno_i(i),3);
        for j = 1:length(barren_list)
          yesterday_x = coordinates_yesterday(barren_list(j),2); yesterday_y = coordinates_yesterday(barren_list(j),3);
          dist_square = (yesterday_x-cellno_x).*(yesterday_x-cellno_x) + (yesterday_y-cellno_y).*(yesterday_y-cellno_y);
          % closest so far
          if dist_square < closest_distance
            closest_distance          = dist_square;
            closest_cellno_yesterday  = barren_list(j);
            closest_cellno_today      = cellno_i(i);
          end
        end
      end
      % update closest
      fprintf(['\n                               ' str3(cellno) ' in fr' str3(yesterdayFrameNum) ' is linked to ' num2str(nr_occur) ' cells in the next frame (' str3(frameNum) '). Corrected: ' num2str(trackLink(closest_cellno_today,:)) ' to ' num2str([closest_cellno_yesterday 0 0 closest_cellno_today])]);
      trackLink(closest_cellno_today,:) = [closest_cellno_yesterday 0 0 closest_cellno_today];

      % recalc
      cellno_i = find(trackLink(:,1)==cellno);
      nr_occur = length(cellno_i);
    end
  end
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % CORRECT FORMAT OF DIVIDING CELLS IN TRACKING
  %------------------------------------------------------------------------
  % now loop over all cellno in yesterday, if any occurs twice, divide
  for cellno = cellnos_yesterday'
    if cellno==0, continue; end % background
    
    % how often does this cellno come back?
    cellno_i = find(trackLink(:,1)==cellno);
    nr_occur = length(cellno_i);
    
    % if cellno occurs 2 times, cell has divided, and trackLink should be adjusted
    if nr_occur==2
      % need to divide
      cellno_today_1 = trackLink(cellno_i(1),4);
      cellno_today_2 = trackLink(cellno_i(2),4);
      trackLink(cellno_i(1),:) = [0 cellno 0 cellno_today_1];
      trackLink(cellno_i(2),:) = [0 0 cellno cellno_today_2];
    end
  end
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % WRITE TRACKING RESULTS TO trackOutputFile
  %------------------------------------------------------------------------
  % clean up any existing output
  if exist(trackOutputFile) == 2
    delete(trackOutputFile)
  end

  % Open trackOutputFile
  fid = fopen(trackOutputFile,'wt');

  % loop over results and print to file
  for i = 1:length(trackLink(:,1))
    fprintf(fid, '%u %u %u %u\n', trackLink(i,1), trackLink(i,2), trackLink(i,3), trackLink(i,4)); 
  end
  
  % Close trackOutputFile
  fclose(fid);
  fprintf(1,' -> finished\n');
%------------------------------------------------------------------------
end  
    
    
%--------------------------------------------------------------------------
% PUT THE TRACKED FRAME PAIRS TOGETHER IN SCHNITZCELLS FILE
%--------------------------------------------------------------------------
MW_makeSchnitzFileFromTracking(p, opts);
%--------------------------------------------------------------------------
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns Lc_fullsize_centered  = cell-centered Lc in fullsize image
%         Lc_fullsize           = Lc in fullsize image
%         Lc                    = original Lc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lc_fullsize_centered, Lc_fullsize, Lc] = loadLcData(segFile);
load(segFile);
Lc_fullsize           = zeros(phaseFullSize);
Lc_fullsize_centered  = zeros(phaseFullSize);

% if segmentation is not approved, get unapproved segmentation 
if ~(exist('Lc') == 1) 
  Lc = LNsub;
end

% put Lc back in original location in fullsize
Lc_fullsize(rect(1):rect(3), rect(2):rect(4)) = Lc;

% get weighted center of cells in Lc_fullsize
[fy, fx]= find(Lc_fullsize>0);
center_cells_x = round( mean(fx) ); % before: round( (max(fx)+min(fx))/2 );
center_cells_y = round( mean(fy) ); % before: round( (max(fy)+min(fy))/2 );

% determine offset
offset_x = round( size(Lc_fullsize,2)/2 - center_cells_x);
offset_y = round( size(Lc_fullsize,1)/2 - center_cells_y);

% center Lc_fullsize into Lc_centered % MW edit 2015/08
% Admininstration I
minfy = min(fy);
maxfy = max(fy);
minfx = min(fx);
maxfx = max(fx);
% Admininstration II
starty  = minfy+offset_y;
endy    = maxfy+offset_y;
startx  = minfx+offset_x;
endx    = maxfx+offset_x;
% Now it could fall outside the image, in that case, just crop it
% This is a bit of a boilerplate solution.. - MW
cropErrorFlag = 0;
if starty<1
    warning('starty<1; Will not perform shift!');
    howmuchoutofplace=1-starty
    starty=1
    %minfy = minfy+howmuchoutofplace
    % Resetting values
    cropErrorFlag = 1;
end
if startx<1
    warning('startx<1; Will not perform shift!');
    howmuchoutofplace=1-startx
    startx=1
    %minfx = minfx+howmuchoutofplace
    cropErrorFlag = 1;
end
if endy>phaseFullSize(1)
    warning('endy>phaseFullSize(1); Will not perform shift!');
    howmuchoutofplace=phaseFullSize(1)-endy
    endy=phaseFullSize(1)
    %maxfy = maxfy+howmuchoutofplace
    cropErrorFlag = 1;
end
if endx>phaseFullSize(2)
    warning('endx>phaseFullSize(1); Will not perform shift!');
    howmuchoutofplace=phaseFullSize(2)-endx
    endx=phaseFullSize(2)
    %maxfx = maxfx+howmuchoutofplace
    cropErrorFlag = 1;
end
% Remark MW 2016/04
% When new center is such that to-select-area falls outside original image
% when performing shift, don't perform the shift..
% Previous solution was to crop, but that resulted into issues.

% Do the centering
if ~cropErrorFlag
    Lc_fullsize_centered( starty:endy, startx:endx ) = Lc_fullsize( minfy:maxfy, minfx:maxfx );
else
    Lc_fullsize_centered=Lc_fullsize;
end;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centroids = getCentroidsOfCell(Lc, cellno);
% extract subcell
cell = +(Lc == cellno);
[fx, fy]= find(cell);
extra= 2;
xmin= max(min(fx) - extra, 1);
xmax= min(max(fx) + extra, size(cell,1));
ymin= max(min(fy) - extra, 1);
ymax= min(max(fy) + extra, size(cell,2));
subcell= cell(xmin:xmax, ymin:ymax);

% find thin
thin = bwmorphmelow(subcell, 'thin', inf);

% reduce thin until only 2 spur points
% get # spur points
spurs = thin & ~bwmorphmelow(thin, 'spur', 1);
[sx,sy] = find(spurs>0);
while length(sx)>2
  % try to remove forked ends
  thin = bwmorphmelow(thin, 'spur', 1);

  % get # spur points
  spurs = thin & ~bwmorphmelow(thin, 'spur', 1);
  [sx,sy] = find(spurs>0);
end

% get coordinates of thin
[px, py]= walkthin(thin);

% find centroids
if length(px)>0
  % centroids is are thirds of thin
  cen(1) = round(length(px)/4);
  cen(2) = round(length(px)/2);
  cen(3) = round(3*length(px)/4);
  centroids(1)   = cellno;
  centroids(2:3) = [(px(cen(2)) + xmin - 1) (py(cen(2)) + ymin - 1)];
  centroids(4:5) = [(px(cen(1)) + xmin - 1) (py(cen(1)) + ymin - 1)];
  centroids(6:7) = [(px(cen(3)) + xmin - 1) (py(cen(3)) + ymin - 1)];
else
  disp('Thin to short?');
  centroids = [cellno round(mean(fx)) round(mean(fy)) round(mean(fx)) round(mean(fy)) round(mean(fx)) round(mean(fy))];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centroids_out = normalizeCentroids(centroids_in, Lc);
[fx, fy]= find(Lc>0);
center_x = (max(fx)+min(fx)) /2 ;
center_y = (max(fy)+min(fy)) /2 ;
height = max(fx)-min(fx); % swapped height and width on 090225 
width = max(fy)-min(fy);

desired_width = size(Lc,2) - 50;
desired_height = size(Lc,1) - 50;

scaling_y = desired_width / width;
scaling_x = desired_height / height;

for i = 1:size(centroids_in,1)
  centroids_out(i,:) = [  centroids_in(i,1) ...
                          center_x + scaling_x * (centroids_in(i,2) - center_x) ...
                          center_y + scaling_y * (centroids_in(i,3) - center_y) ...
                          center_x + scaling_x * (centroids_in(i,4) - center_x) ...
                          center_y + scaling_y * (centroids_in(i,5) - center_y) ...
                          center_x + scaling_x * (centroids_in(i,6) - center_x) ...
                          center_y + scaling_y * (centroids_in(i,7) - center_y) ];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Returns Lc          = original Lc (size differs for each frame)
% %         Lc_centered = Lc reset in original full size image
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Lc_centered, Lc] = loadCenteredSegData(segFile);
% load(segFile);
% Lc_centered = zeros(phaseFullSize);
% 
% % if segmentation is not approved, get unapproved segmentation 
% if ~(exist('Lc') == 1) 
%   Lc = LNsub;
% end
% 
% % put back in center 
% off_set_x = uint16( (0.5+phaseFullSize(1)-rect(1)-rect(3))/2 );
% off_set_y = uint16( (0.5+phaseFullSize(2)-rect(2)-rect(4))/2 );
% 
% % put Lc back in original location
% Lc_centered(rect(1)+off_set_x:rect(3)+off_set_x, rect(2)+off_set_y:rect(4)+off_set_y) = Lc;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Lc_centered, Lc] = loadCenteredSegData2(segFile);
% % get data
% [Lc_fullsize, Lc] = loadCenteredSegData(segFile);
% 
% % get center of Lc_fullsize
% [fy, fx]= find(Lc_fullsize>0);
% center_Lc_x = round( (max(fx)+min(fx))/2 );
% center_Lc_y = round( (max(fy)+min(fy))/2 );
% 
% % determine offset
% offset_x = round( size(Lc_fullsize,2)/2 - center_Lc_x);
% offset_y = round( size(Lc_fullsize,1)/2 - center_Lc_y);
% 
% % center Lc_fullsize into Lc_centered
% Lc_centered = zeros(size(Lc_fullsize));
% Lc_centered( min(fy)+offset_y:max(fy)+offset_y, min(fx)+offset_x:max(fx)+offset_x ) = Lc_fullsize( min(fy):max(fy), min(fx):max(fx) );
% 
% % DEBUGGING : show change
% if false
%   disp(['x1: border left center right border : 0 ' str3(min(fx)) ' ' str3(center_Lc_x) ' ' str3(max(fx)) ' ' str3(size(Lc_fullsize,2))]);
%   disp(['y1: border left center right border : 0 ' str3(min(fy)) ' ' str3(center_Lc_y) ' ' str3(max(fy)) ' ' str3(size(Lc_fullsize,1))]);
% 
%   [fy2, fx2]= find(Lc_centered>0);
%   center_Lc_x2 = round( (max(fx2)+min(fx2))/2 );
%   center_Lc_y2 = round( (max(fy2)+min(fy2))/2 );
% 
%   disp(['x2: border left center right border : 0 ' str3(min(fx2)) ' ' str3(center_Lc_x2) ' ' str3(max(fx2)) ' ' str3(size(Lc_centered,2))]);
%   disp(['y2: border left center right border : 0 ' str3(min(fy2)) ' ' str3(center_Lc_y2) ' ' str3(max(fy2)) ' ' str3(size(Lc_centered,1))]);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
