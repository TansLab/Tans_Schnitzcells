function p = NW_tracker_centroid_vs_area(p, varargin)
% NW_tracker_centroid_vs_area tracks cells across frames in one complete movie.
%
% 1) Each cell is tracked between frame pairs by minimizing the distance of
% its current backbone centroid (center of 'thin' line) to all possible
% areas of cells in the previous image.
% 2) Barren cells and cells with more than 2 offspring are automatically
% accounted for.
% 3) Tracked cells which display suspiciously large increase/decrease in area
% are reconnected by evaluating a cost function
%    ---- optimal parameters for this cost function may alter -----
%    ---- and can be adjusted below in this main function!!   -----
%
% The function scaffold and data struture is copied from DJK_tracker_djk, 
% however the tracking method and tracking-checks are renewed.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT:
%   p
%
% OPTIONAL INPUT:
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
%                 previously performed for those frames (that is, the
%                 segmentation file is older than the tracking file) 
%   SwapcounterMax   How many subsequent rounds shall possibly mistracked
%                 cells be pairwise swapped. default: 1
%                 if swapped cells remain, increase this number and run
%                 tracker again over this frame pair ('overwrite')
%
% NOTE:
% For tracking correction (cell swapping, division association) several
% parameters are hard coded and can be changed within the script.
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%  schnitzcells   schnitzcells is a 1xNS struct array (NS = number of "schnitz" 
%                 tracks).  This struct array describes the lineage of all 
%                 cells in a movie.  A "schnitz" is a track that begins at 
%                 the first frame a distinct cell appears in, and ends at the 
%                 last frame that distinct cell appears in.  When a cell 
%                 divides its track ends, and two new tracks being at the 
%                 first frame its two children appear in.  A schnitz (track) 
%                 is numbered by its position (index) within the struct array. 
%                 A minimal schnitz structure is described by P (the schnitz 
%                 number of it's parent), children schnitz numbers D and E,
%                 sister schnitz number S, frames (an array) listing each frame 
%                 this schnitz occurs in, cellno (an array) listing each 
%                 cell number of the tracked cell within each of the frames, 
%                 and N, the number of frames that distinct cell appears.
% frame1-frame2-trackfiles, e.g. pos5crop-djk-output-467-to-468.txt: Links of
%                 all tracked cell numbers between frame1 & frame2.
%                 Structure:
%                 cellno1  0         0        cellno2
%                 0        parent1   0        daughter2_1
%                 0        0         parent1  daughter2_2
%                 
%                 
%--------------------------------------------------------------------------

% CONSTANTS THAT ARE SET IN THE FUNCTION BELOW-----------------------------
% ******** ADJUST PARAMETERS ************************************************ 
CUTOFFAREA=200; % which area change is suspicious?
CUTOFFDISPLACEMENT=15; % which displacement (of non-dividing cells) is suspicious?
% weighing parameters for cost functions:
% -- if e.g. normalAreaChange=100 and
% normalDisplacement(centroid<->area)=20, then a total area increase of
% 100px at a track link has the same penalty as a total displacement of 20px.
NORMALAREACHANGEDIV=100; %(px) % 2x cell area growth
NORMALDISPLACEMENTDIV=20; %(px) % 3x displacement of centroid<->area !
NORMALAREACHANGENONDIV=100; %(px)
NORMALDISPLACEMENTNONDIV=20; %(px) 
% ***************************************************************************
%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1;
functionName = 'NW_tracker_centroid_vs_area';

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
% if # of swapping cell attempts is not given, try only once
if ~existfield(p,'SwapcounterMax')
  p.SwapcounterMax=1;
end

% Keep only the frames in the range that contain a corrected segmentation (unless we're tracking un-checked frames)
manualRangeChecked = [];
for frameNum = p.manualRange
  clear Lc 
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

%**********************************************************************************************************
%--------------------------------------------------------------------------
% CARRY OUT TRACKING -> LOOP OVER FRAME PAIRS
%--------------------------------------------------------------------------
%**********************************************************************************************************frameNum = 0;
%dataLeft_previousRound = false; not applicable any more since areas <-> centroids
for count = 2:length(p.manualRange);
    clear allcentroidstoday allcentroidsyesterday;
    
  %frameNum_previousRound = frameNum; not applicable any more: area<->centroids
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
  clear Lc_today_fullsize Lc_yesterday_fullsize % unnecessary
  % rescale yesterday image to the size of the today image
  Lc_yesterday_fullsize_resc_centered=rescaleLcData(Lc_yesterday,Lc_today, ...
      Lc_yesterday_fullsize_centered); %use subimages as input
  % CAREFUL: Use the following images to obtain coordinates etc:
            % Lc_today_fullsize_centered     and
            % Lc_yesterday_fullsize_RESC_centered   !!
  
  
  % - - - - - - - - - - - - - - - - - - -
  % SAVE SOME PORPERTIES FOR LATER ANALYSIS STEPS (DJK_analyze_tracking) . 
  % Not directly needed for tracking any more (was needed in DJK_tracker_djk). 
  % Data is added to schnitzes in recalc_schnitz below.
  % - - - - - - - - - - - - - - - - - - -
  rp  = regionprops(Lc_today,'Centroid','Orientation','MajorAxisLength','Area');
  rp2 = regionprops(Lc_today_fullsize_centered,'Centroid');
  num_pts = size(rp,1);
  for j=1:num_pts
    Points(j).cenx      = rp(j).Centroid(1);
    Points(j).ceny      = rp(j).Centroid(2);
    Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
    Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
    Points(j).ang       = rp(j).Orientation;
    Points(j).len       = rp(j).MajorAxisLength;
    Points(j).areapx       = rp(j).Area;  %NW 2013-12
    Points(j).cellno    = j;
  end
  opts{count}=Points(1:num_pts);
  
  % in case of first run, one extra
  if count == 2
    rp  = regionprops(Lc_yesterday,'Centroid','Orientation','MajorAxisLength','Area');
    rp2 = regionprops(Lc_yesterday_fullsize_centered,'Centroid'); % don't use rescale!
    num_pts = size(rp,1);
    for j=1:num_pts
      Points(j).cenx    = rp(j).Centroid(1);
      Points(j).ceny    = rp(j).Centroid(2);
      Points(j).cenx_cent = rp2(j).Centroid(1); % DJK 090410
      Points(j).ceny_cent = rp2(j).Centroid(2); % DJK 090410
      Points(j).ang     = rp(j).Orientation;
      Points(j).len     = rp(j).MajorAxisLength;
      Points(j).areapx       = rp(j).Area;  %NW 2013-12
      Points(j).cellno  = j;
    end
    opts{count-1}=Points(1:num_pts);
  end
  % - - - - - - - - - - - - - - - - - - -
            
  % tell what is happening. continue if frames are already tracked.
  fprintf([' * frame pair ',str3(yesterdayFrameNum),' -> ', str3(frameNum) ' : ']);
  if ~needToTrack
    %dataLeft_previousRound = false;  %not applicable any more (area<->centroid)
    fprintf(1,' -> Skipping, cause seg older than previous tracking (use p.overwrite=1 to redo)\n');
    continue
  end          
            
  % - - - - - - - - - - - - - - - - - - -
  % IF FRAMES STILL NEED TO BE TRACKED:
  % - - - - - - - - - - - - - - - - - - -
  
  %cell numbers in Lc frames (ASSURE THAT CELLNO's ARE 1:max2(Lc), if gaps in cell numbers: return error )
  all_cellno_today=unique(Lc_today_fullsize_centered)'; %transpose!
  if all_cellno_today(1)==0 % remove background (=0)
      all_cellno_today=all_cellno_today(2:end);
  end
  all_cellno_yesterday=unique(Lc_yesterday_fullsize_resc_centered)'; %transpose!
  if all_cellno_yesterday(1)==0 % remove background (=0)
      all_cellno_yesterday=all_cellno_yesterday(2:end);
  end
  if ~isequal(all_cellno_today,1:max2(Lc_today_fullsize_centered))
      errormessage=['Problem with cell numbers! Segmentation file ' num2str(frameNum) ...
          ' does not contain cell numbers 1:n. Tracker can be adjusted to deal with jumps ' ...
           ' in cell numbers. However, a segmentation issue seems to be present.'];
      error(errormessage);
  end
  if ~isequal(all_cellno_yesterday,1:max2(Lc_yesterday_fullsize_resc_centered))
      errormessage=['Problem with cell numbers! Segmentation file ' num2str(yesterdayFrameNum) ...
          ' does not contain cell numbers 1:n. Tracker can be adjusted to deal with jumps ' ...
           ' in cell numbers. However, a segmentation issue seems to be present.'];
      error(errormessage);
  end
  
  %initialize vector for storage of tracking connections
  trackLink=zeros(0,4);
            
  % - - - - - - - - - - - - - - - - - - -
  % GET AREA INFO OF YESTERDAY CELLS AND BACKBONE('THIN')-CENTROIDS OF TODAY CELLS
  % - - - - - - - - - - - - - - - - - - -
  % get coordinates of all pixels of cell areas in yesterday image
  allpixellistYesterday=regionprops(Lc_yesterday_fullsize_resc_centered,'PixelList'); 
  %Pixellist(:,1): Columns
  %Pixellist(:,2): Rows
  
  % get backbone centroids of today cells
  allcentroidstoday=zeros(max(all_cellno_today),2);
  for cellnorun=all_cellno_today
      cellnoruncentroid=getCentroidsOfCell(Lc_today_fullsize_centered, cellnorun);
      % only 1/2 centroid (and not 1/4 and 3/4) is needed
      cellnoruncentroid=cellnoruncentroid([2 3]); %(y,x)=(row,column)
      allcentroidstoday(cellnorun,:)=cellnoruncentroid;
  end
  DEBUGFINDCENTROID=0;
  if DEBUGFINDCENTROID
      centrim=zeros(size(Lc_today_fullsize_centered));
      for i=all_cellno_today
        centrim(allcentroidstoday(i,1),allcentroidstoday(i,2))=1;
      end
      centrim=imdilate(centrim,strel('disk',5));
      combinedim=(Lc_today_fullsize_centered>0)+2*centrim;
      figure
      imagesc(combinedim)
      title('today image with backbone centroids.')
  end
  % - - - - - - - - - - - - - - - - - - -
  %------------------------------------------------------------------------
  
  
  %------------------------------------------------------------------------
  % LOOP OVER ALL CELLS AND TRACK
  %------------------------------------------------------------------------
  
  for cellno_today=all_cellno_today
      currentcentroid=allcentroidstoday(cellno_today,:); %(row, column)
      centroidRow=currentcentroid(1);
      centroidCol=currentcentroid(2); 
      % create storage vector for all distances to yesterday cells
      mindistVec=zeros(length(all_cellno_yesterday),1);
      
      %calculate distance to each yesterday cell
     for cellno_yesterday=all_cellno_yesterday
         AreaPixelListRow=allpixellistYesterday(cellno_yesterday).PixelList(:,2);
         AreaPixelListCol=allpixellistYesterday(cellno_yesterday).PixelList(:,1);
         Celldist=GetDistance_Point_from_Area(AreaPixelListRow, AreaPixelListCol, centroidRow,centroidCol);
         mindistVec(cellno_yesterday)=Celldist;
         clear distfromCellYesterday Celldist
      end
      sortedCellDist=sort(mindistVec);
      tracked_cellno_yesterday=find(mindistVec==sortedCellDist(1));
      % if >=2 yesterday cells have exactly the same distance (can happen due
      % to discrete pixel values choose the first yest-cell which is not
      % linked yet (This selection criterion could still be improved NW2013-12)
      if length(tracked_cellno_yesterday)>1
          tmpcells=tracked_cellno_yesterday;
          tracked_cellno_yesterday=[];
          for cell_yest=tmpcells'
              % yest-cell is not connected yet
              if isempty(find(trackLink(:,1),cell_yest))
                  tracked_cellno_yesterday=cell_yest;
                  continue
              end
          end
          % no not connected yest-cell found
          if isempty(tracked_cellno_yesterday)
              tracked_cellno_yesterday=tmpcells(1);
          end
      end 
      trackLink =[trackLink; tracked_cellno_yesterday 0 0 cellno_today];
      
      clear mindistVec cellidx sortedCellDist tracked_cellno_yesterday cellno_yesterday
  end
  clear cellno_today
  % -----------------------------------------------------------------------------
  
  
  
  
  
% *****************************************************************************************

  %------------------------------------------------------------------------
  % CORRECT ERRORS IN TRACKING: BARREN CELLS
  %------------------------------------------------------------------------
  % loop over all cellno in yesterday, if any does not occur in the
  % today image, correct (create a new link)
  barrencells=[];
  for cellno_yesterday = all_cellno_yesterday
    % how often is the yesterday cell (cell_no_yesterday) connected?
    cellno_children = find(trackLink(:,1)==cellno_yesterday);
    nr_occur = length(cellno_children);
    if nr_occur==0
        barrencells=[barrencells , cellno_yesterday];
    end
  end
  
  if ~isempty(barrencells)
  % IF BARREN CELLS EXIST, WE HAVE TO CORECT FOR IT:
  % 1) get area info of today(!) cells and backbone centroids of yesterday(!)
  % 2) find dividing cells (possible new connections)
  % 3) find an inverse tracking connection similar to the standard tracking
  
  % 1) ** get today areas and yesterday backbone centroids **
  allpixellistToday=regionprops(Lc_today_fullsize_centered,'PixelList');
  allcentroidsyesterday=zeros(max(all_cellno_yesterday),2);
  for cellnorun=all_cellno_yesterday
      cellnoruncentroid=getCentroidsOfCell(Lc_yesterday_fullsize_centered, cellnorun);
      % only 1/2 centroid (and not 1/4 and 3/4) is needed
      cellnoruncentroid=cellnoruncentroid([2 3]); %(y,x)=(row,column)
      allcentroidsyesterday(cellnorun,:)=cellnoruncentroid;
  end
  
  % 2) ** find dividing cells **
  % loop over cells in yesterday that are linked to >= 2 cells
  today_dividing_list = []; % all daughter cells in today's image which are newly divided and could be connected to the barren cell.
  for cellno_maybeparent = all_cellno_yesterday
       % get the cells they link to
       linked_cells = find(trackLink(:,1)==cellno_maybeparent);
       if length(linked_cells) >= 2 % at least two 'today cell' are linked to yesterday cell 'cellno_maybeparent'
                                    % -> cell division
           % now add the cells they link to today_dividing_list
           today_dividing_list = [today_dividing_list linked_cells']; 
       end
  end
      
  % 3) ** connect barren cells **
  for cellno_yesterday = barrencells
      if length(today_dividing_list) > 0 % -> some divisions occured
              % calc distances between yesterday (barren) centroids and todays areas (daughter cells) and pick closest
              currentcentroidYest=allcentroidsyesterday(cellno_yesterday,:); %(row, column) Careful: Yesterday, not Today in this case!!
              centroidRowYest=currentcentroidYest(1);
              centroidColYest=currentcentroidYest(2); 
              
              % create storage vector for distance to yesterday cells
              mindistVec=[];

              % loop over all possible today-cells
              for cellno_today=today_dividing_list % !! Different than the other loops: Here the loop
                                    % runs over today-cells and not yesterday-cells !!
                  % get distance of this today cell with the barren yesterday cell:
                  AreaPixelListRowToday=allpixellistToday(cellno_today).PixelList(:,2);
                  AreaPixelListColToday=allpixellistToday(cellno_today).PixelList(:,1);
                  Celldist=GetDistance_Point_from_Area(AreaPixelListRowToday, AreaPixelListColToday, ... 
                      centroidRowYest,centroidColYest);
                  mindistVec=[mindistVec; Celldist];
                  clear Celldist 
              end
              sortedCellDist=sort(mindistVec);
              idx_cellno_today=find(mindistVec==sortedCellDist(1));
              % if same minimal distance to more than one cell: use the one
              % with the lowest cell number (criterion could be improved
              % NW2013-12)
              if length(idx_cellno_today)>1
                  idx_cellno_today=idx_cellno_today(1)
              end
              tracked_cellno_today=today_dividing_list(idx_cellno_today);
              %check if new yesterday cell is still connected to at least
              %two today cells. otherwise remove from today_dividing_list
              % NW 2014-01
              remove_cells_from_divlist=[]; % cells which are not result of division any more should be removed from dividing list
              old_parent_idx=find(trackLink(:,4)==tracked_cellno_today); %should be same as tracked_cellno_today
              old_parent=trackLink(old_parent_idx);
              daughters=find(trackLink(:,1)==old_parent); % how many daughters were there?
                               
              
              clear mindistVec cellidx sortedCellDist;

              % update barren cell with new connection to closest (daughter) cell
              fprintf(['\n                ' str3(cellno_yesterday) ' in fr' str3(yesterdayFrameNum) ... 
                    ' was not linked. Corrected: ' num2str(trackLink(tracked_cellno_today,:)) ' to ' ...
                    num2str([cellno_yesterday 0 0 tracked_cellno_today])]);
              % recalc
              trackLink(tracked_cellno_today,:) = [cellno_yesterday 0 0 tracked_cellno_today];
              
              % if only 2 daughters before rewiring, just one track is left
              % after rewiring -> remove from today_dividing list
              if length(daughters)<=2
                  remove_cells_from_divlist=trackLink(daughters,4);
                  fprintf([' \n             Will remove today cells ' num2str(daughters') ' from dividing list.'])
              % if 3 links existed, remove only the daughter which is now
              % newly connected to the former barren cell, i.e. is part of
              % a non-div connection
              else
                  remove_cells_from_divlist=trackLink(tracked_cellno_today,4); %NW 2014-05
                                            % trackLink(tracked_cellno_today,4)
                                            % should be the same as
                                            % 'tracked_cellno_today'
                  fprintf([' \n             Will remove today cells ' num2str(trackLink(tracked_cellno_today,4)) ' from dividing list.'])
              end

              % MW 2015/06
              % now remove cells which are not dividing any more from the list:
              today_dividing_list=today_dividing_list(~ismember(today_dividing_list,remove_cells_from_divlist));
             
        else
            fprintf(['\n                  ' str3(cellno_yesterday) ' in fr' ... 
                str3(yesterdayFrameNum) ' was not linked. Could not be corrected, cause no dividing cells!']);
      end

      
    end
  end
  %------------------------------------------------------------------------
  
  
  
  %------------------------------------------------------------------------
  % CORRECT ERRORS IN TRACKING: MULTIPLE LINKS
  %------------------------------------------------------------------------
  % now loop over all cellno in yesterday, if any is connected > twice
  % (that is e.g. three daughters), correct
  for cellno_yesterday = all_cellno_yesterday
    % how often is the yesterday cell (cell_no) connected?
    cellno_children = find(trackLink(:,1)==cellno_yesterday); % a list of all today-cell-numbers (daughters) 
                                                    % which are linked to the same yesterday cell
    nr_occur = length(cellno_children); % number of daughters    
    % if cellno occurs more than 2 times (>2 daughters), it must be linked to extra cell(s) 
    % -> correct for this
    if nr_occur>2 
      fprintf(['\n                               ' str3(cellno_yesterday) ' in fr' str3(yesterdayFrameNum) ' is linked to ' num2str(nr_occur) ' cells in the next frame (' str3(frameNum) ')']);
    end
    
    % do correction, until nr_occur==2
    while nr_occur>2 
      % find cells in yesterday which can be used for new connection to the
      % surplus daughter cells:
      % get list of cells in yesterday frame which are linked to no cells
      barren_list = []; % barren_list: all cellnumbers of yesterday which are linked to at
                        % most one cell of today
      for cellno_yes = all_cellno_yesterday
        if length(find(trackLink(:,1)==cellno_yes)) == 0
          barren_list = [barren_list cellno_yes];
        end
      end
      % if no more barren cells, will continue with cells linked to 1 other cell
      if length(barren_list) == 0
        for cellno_yes = all_cellno_yesterday
          if length(find(trackLink(:,1)==cellno_yes)) == 1
            barren_list = [barren_list cellno_yes];
          end
        end
      end
      
      % calc distances between surplus today cells and suitable yesterday cells and pick closest
      closest_distance          = realmax('double'); % of all (typically 3) daughter cells
                                                    % with all possible yesterday cells
      closest_cellno_yesterday  = 1;
      closest_cellno_today      = 1;
      for i = 1:nr_occur % loop over all (typically =3) daughter cells
          % here: get centroids of todays daughter cells
         currentcentroid=allcentroidstoday(cellno_children(i),:); %(row, column)
         centroidRow=currentcentroid(1);
         centroidCol=currentcentroid(2); 
        for j = 1:length(barren_list) %loop over possible yesterday cells and get distance from yest areas 
            AreaPixelListRow=allpixellistYesterday(barren_list(j)).PixelList(:,2);
            AreaPixelListCol=allpixellistYesterday(barren_list(j)).PixelList(:,1);
            Celldist=GetDistance_Point_from_Area(AreaPixelListRow, AreaPixelListCol, centroidRow,centroidCol);
            if Celldist < closest_distance
                closest_distance          = Celldist; %choses closest distance among all (typically 3) daughters
                closest_cellno_yesterday  = barren_list(j); %yesterday cell corresponding to the closest distance
                closest_cellno_today      = cellno_children(i); %today (daughter) cell which will be associated to a new yest-cell
            end
        end
      end
      
      % update closest daughter cell with new yesterday cell      
      fprintf(['\n                               ' str3(cellno_yesterday) ' in fr' str3(yesterdayFrameNum) ... 
          ' is linked to ' num2str(nr_occur) ' cells in the next frame (' str3(frameNum) '). Corrected: ' ... 
          num2str(trackLink(closest_cellno_today,:)) ' to ' num2str([closest_cellno_yesterday 0 0 closest_cellno_today])]);
      trackLink(closest_cellno_today,:) = [closest_cellno_yesterday 0 0 closest_cellno_today];
      % recalc trackLink
      cellno_children = find(trackLink(:,1)==cellno_yesterday);
      nr_occur = length(cellno_children);
    end
  end
  %------------------------------------------------------------------------
  
  
  % -----------------------------------------------------------------------
% THE FOLLOWING ERROR CHECKS TAKE CARE OF TRACKING ERRORS TYPICALLY
% ENCOUNTERED IN THE CURRENT EXPERIMENTS -> IF ERROR CLASSES IN THE EXPERIMENT
% OF YOUR CHOICE ARE NOT ACCOUNTED FOR, A SPECIFIC CORRECTION CAN BE INCLUDED
% HERE (e.g. moving groups of cells (>2) might need an extra correction)
% NW 2013-12
% -----------------------------------------------------------------------


% ------------------------------------------------------------------------
% CORRECT SUSPICIOUS JUMPS IN CELL AREAS:
% EITHER WRONGLY TRACKED DIVISIONS OR CELL SWAPPINGS
% NB: CHECK FOR SWAPS COULD BE IMPROVED BY INTRODUCING DISPLACEMENT AS
% SECOND CRITERION (APART FROM AREA CHANGE)
% ------------------------------------------------------------------------
% all threshold values in these subfunctions can be modified below.
% steps: 
% 1) Find dividing cells
% 2) Get areas of all (dividing and non-dividing) cells
% 3) Calculate Delta(Area) for all cells (dividing cells: account for 2
%    daughters) and select strange changes (e.g.: more than 150px)
% XY) find dividing cells with suspicuously high displacement
% 4) Seperate strange tracks into dividing and non-dividing cells
% 5) Try relinking the dividing cells with another non-div cell
%    Use a costfunction that accounts for area change & displacement as
%    tracking evaluation
% 6) Find suspicuous non-dividing cells (area change or large movement).
%     Try relinking the non-dividing cells


% *****************************************************************************************
% *************** STRANGE DIVISIONS **************
% *****************************************************************************************

%Brute force:
% Run whole swapping x times so that triplett swaps etc can be detected
% Should be improved by considering a matrix of possible cell swaps, not
% only pairs.

 for Swapcounter=1:p.SwapcounterMax

% 1) find yesterday cells which are linked to 2 daugthers (>2 already
% corrected for above) -> cell divisions
parent_vec=[];  
for cellno_maybeparent = all_cellno_yesterday
       % get the cells they link to
       linked_cells = find(trackLink(:,1)==cellno_maybeparent);
       if length(linked_cells) == 2 % two 'today cell' are linked to yesterday cell 'cellno_maybeparent'
                                    % -> cell division
           parent_vec=[parent_vec; cellno_maybeparent];  
       end
end

% -----------------------------------------------------
% Find strange area changes for non-dividing cells
% -----------------------------------------------------

% 2) get areas of all cells (today & yesterday)
allareasyesterday=regionprops(Lc_yesterday_fullsize_resc_centered,'Area'); %use rescaled (?)
allareasyesterday=[allareasyesterday.Area];
allareastoday=regionprops(Lc_today_fullsize_centered,'Area');
allareastoday=[allareastoday.Area];

% 3) calculate area increase of all tracked cells (if division:
% accounts for 2 daughters)
allareas_trackedtotoday=zeros(size(allareasyesterday));  % idx -> cell number of yesterday
for cellno_yest=all_cellno_yesterday
    linked_cells=find(trackLink(:,1)==cellno_yest);
    allareas_trackedtotoday(cellno_yest)=sum(allareastoday(linked_cells));
    clear linked_cells
end
clear cellno_yest
deltaArea=allareas_trackedtotoday-allareasyesterday;
%find yesterday cells with strange change in area (e.g.: >150px,<-150px)
yestcellsStrangeArea=find(abs(deltaArea)>=CUTOFFAREA); %-> yest cell nrs with strange area change

% -----------------------------------------------------
% Find Strangely Moving Cells (Displacement (centroid-centroid [not area]))
% -----------------------------------------------------
% 1) if yesterday centroid vector doesn't exist, create it
if ~exist('allcentroidsyesterday')
    allcentroidsyesterday=zeros(max(all_cellno_yesterday),2);
      for cellnorun=all_cellno_yesterday
          cellnoruncentroid=getCentroidsOfCell(Lc_yesterday_fullsize_centered, cellnorun);
          % only 1/2 centroid (and not 1/4 and 3/4) is needed
          cellnoruncentroid=cellnoruncentroid([2 3]); %(y,x)=(row,column)
          allcentroidsyesterday(cellnorun,:)=cellnoruncentroid;
      end
end
% 2) calculate displacements of centroids
allcentroids_trackedtotoday=zeros(size(allcentroidsyesterday));  % idx -> cell number of yesterday
for cellno_yest=all_cellno_yesterday
    linked_cells=find(trackLink(:,1)==cellno_yest);
    % if the cell has dividided take the mean coordinates of all daughter
    % centroids
    if length(linked_cells)==1 %not divided
        allcentroids_trackedtotoday(cellno_yest,:)=allcentroidstoday(linked_cells,:);
    elseif length(linked_cells)>1 % divided
        allcentroids_trackedtotoday(cellno_yest,:)=mean(allcentroidstoday(linked_cells,:));
    end % if no link at all: entry is [0 0]
end
clear cellno_yest
dummy=(allcentroids_trackedtotoday-allcentroidsyesterday);
deltaDisplacement=sqrt(dummy(:,1).^2+dummy(:,2).^2);
%find yesterday cells with strange movement (e.g.: >30px)
% still includes div & non-div cells
yestcellsStrangeDisplacement=find(deltaDisplacement>=CUTOFFDISPLACEMENT); %-> yest cell nrs with strange area change

% --------------------------------------------------
% ------- combine strange yest. cells and separate into div/non-div cells ---
% --------------------------------------------------
yestcellsStrange_Area_or_Displacement=unique([yestcellsStrangeArea, yestcellsStrangeDisplacement']);

if ~isempty(yestcellsStrange_Area_or_Displacement)
    % 4) find dividing & non-dividing yesterday cells which have strange change in area
    idxStrangeDiv=find(ismember(yestcellsStrange_Area_or_Displacement,parent_vec)); 
    yestcellsStrange_Area_or_Displacement_Div=yestcellsStrange_Area_or_Displacement(idxStrangeDiv);
    idxStrangeNonDiv=find(~ismember(yestcellsStrange_Area_or_Displacement,parent_vec));
    yestcellsStrange_Area_or_Displacement_NonDiv=yestcellsStrange_Area_or_Displacement(idxStrangeNonDiv);
end


% ------------------------------------------------------
% Attempt correction of wrong division events
% ------------------------------------------------------

% div & non-div cells exist
if ~isempty(yestcellsStrange_Area_or_Displacement)
    if ~isempty(yestcellsStrange_Area_or_Displacement_Div) & ~isempty(yestcellsStrange_Area_or_Displacement_NonDiv)
    
    % 5) if strange divisions occured, correct first tracking in divisions
    
        % create cost function for possible tracking links. THIS IS VERY
        % EMPIRICAL AND MIGHT HAVE TO BE ADJUSTED!!
        % also no renormalization of areas & displacements
        % normalization parameters: defined above
        %blubb mycostfunction = @ (deltaarea1, deltaarea2, deltadist1, deltadist2, deltadist3) ...
            %(deltaarea1 + deltaarea2)/normalAreaChangeDiv + (deltadist1+deltadist2+deltadist3)/normalDisplacementDiv;
            
            % **** IMPORTANT ****
            % Division event: deltaarea2, deltadist2, deltadist3.
            % KEEP CORRECT ORDER OF INPUT ARGUMENTS!
            % Divions event is weight stronger
            % ******************
        mycostfunction = @ (deltaarea1, deltaarea2, deltadist1, deltadist2, deltadist3) ...
            (deltaarea1 + deltaarea2)/NORMALAREACHANGEDIV + (deltadist1+deltadist2+deltadist3)/NORMALDISPLACEMENTDIV+ ...
            + deltaarea2/NORMALAREACHANGENONDIV *10 + (deltadist2+deltadist3)/NORMALDISPLACEMENTNONDIV *10; %overemphasize divsion event
        
        % input: all area changes (2), all displacements (3)
        % output: cost of tracking connection

        % loop over all strange-area dividing cells
            for y1=yestcellsStrange_Area_or_Displacement_Div
                CellsCorrected=0;
                cellnos_todayDiv=find(trackLink(:,1)==y1);
                t1=cellnos_todayDiv(1); t2=cellnos_todayDiv(2); % tracked division cells
                % loop over all strange-area non-dividing cells (possible
                % switch partners)
                for y2=yestcellsStrange_Area_or_Displacement_NonDiv
                    % get all 3 cells of today which are linked to (cellno_yestDiv and cellno_yestNonDiv)
                    cellnos_todayNonDiv=find(trackLink(:,1)==y2);
                    if length(cellnos_todayNonDiv)~=1 
                        fprintf(' \n Warning: Error in tracking (area) correction of dividing cells.')
                        fprintf('\n cellnos_today_NonDiv should have one link. Will skip. Barren cell?')
                        continue
                    end

                    t3=cellnos_todayNonDiv(1); % tracked non-div cell.
                    % permute over tracking possibilities and find all abs. area
                    % changes & centroid displacements (could be improved with
                    % a loop/matrix)
                    % for convenience and readability, precalculate all
                    % possible abs area changes and centroid<->area
                    % displacements
                    area_y1t1=abs(allareastoday(t1)-allareasyesterday(y1));   % area y1->t1 (no division)
                    area_y1t2=abs(allareastoday(t2)-allareasyesterday(y1));
                    area_y1t3=abs(allareastoday(t3)-allareasyesterday(y1));
                    area_y2t1=abs(allareastoday(t1)-allareasyesterday(y2));
                    area_y2t2=abs(allareastoday(t2)-allareasyesterday(y2));
                    area_y2t3=abs(allareastoday(t3)-allareasyesterday(y2));
                    area_y1t1t2=abs(allareastoday(t1)+allareastoday(t2)-allareasyesterday(y1)); % area y1->t1+t2 (division)
                    area_y1t1t3=abs(allareastoday(t1)+allareastoday(t3)-allareasyesterday(y1));
                    area_y1t2t3=abs(allareastoday(t2)+allareastoday(t3)-allareasyesterday(y1));
                    area_y2t1t2=abs(allareastoday(t1)+allareastoday(t2)-allareasyesterday(y2));
                    area_y2t1t3=abs(allareastoday(t1)+allareastoday(t3)-allareasyesterday(y2));
                    area_y2t2t3=abs(allareastoday(t2)+allareastoday(t3)-allareasyesterday(y2));

                    dist_y1t1=GetDistance_Point_from_Area(allpixellistYesterday(y1).PixelList(:,2), allpixellistYesterday(y1).PixelList(:,1), allcentroidstoday(t1,1), allcentroidstoday(t1,2));
                    dist_y1t2=GetDistance_Point_from_Area(allpixellistYesterday(y1).PixelList(:,2), allpixellistYesterday(y1).PixelList(:,1), allcentroidstoday(t2,1), allcentroidstoday(t2,2));
                    dist_y1t3=GetDistance_Point_from_Area(allpixellistYesterday(y1).PixelList(:,2), allpixellistYesterday(y1).PixelList(:,1), allcentroidstoday(t3,1), allcentroidstoday(t3,2));
                    dist_y2t1=GetDistance_Point_from_Area(allpixellistYesterday(y2).PixelList(:,2), allpixellistYesterday(y2).PixelList(:,1), allcentroidstoday(t1,1), allcentroidstoday(t1,2));
                    dist_y2t2=GetDistance_Point_from_Area(allpixellistYesterday(y2).PixelList(:,2), allpixellistYesterday(y2).PixelList(:,1), allcentroidstoday(t2,1), allcentroidstoday(t2,2));
                    dist_y2t3=GetDistance_Point_from_Area(allpixellistYesterday(y2).PixelList(:,2), allpixellistYesterday(y2).PixelList(:,1), allcentroidstoday(t3,1), allcentroidstoday(t3,2));

                    % calculate cost functions of tracking & create an array to easily
                    % recover the cell numbers
                    matrix_cellno_cost=zeros(6,3); % per row: non-div yest cell  --- non-div today cell --- cost
                                                   % 6 tracking options (1x div, 1x non-div)
                    % ********* ALWAYS WRITE THE DIVISION AREA &
                    % DISPLACEMEENT AS 2nd resp 2nd & 3rd *****
                                                   
                    cost_y1t1_y2t2t3=mycostfunction(area_y1t1,area_y2t2t3,dist_y1t1,dist_y2t2,dist_y2t3); % y1->t1, y2->t2&t3
                    matrix_cellno_cost(1,:)=([y1,t1,cost_y1t1_y2t2t3]);
                    cost_y1t2_y2t1t3=mycostfunction(area_y1t2,area_y2t1t3,dist_y1t2,dist_y2t1,dist_y2t3);
                    matrix_cellno_cost(2,:)=([y1,t2,cost_y1t2_y2t1t3]);
                    cost_y1t3_y2t1t2=mycostfunction(area_y1t3,area_y2t1t2,dist_y1t3,dist_y2t1,dist_y2t2);
                    matrix_cellno_cost(3,:)=([y1,t3,cost_y1t3_y2t1t2]);
                    cost_y2t1_y1t2t3=mycostfunction(area_y2t1,area_y1t2t3,dist_y2t1,dist_y1t2,dist_y1t3); % y2->t1, y1->t2&t3
                    matrix_cellno_cost(4,:)=([y2,t1,cost_y2t1_y1t2t3]);
                    cost_y2t2_y1t1t3=mycostfunction(area_y2t2,area_y1t1t3,dist_y2t2,dist_y1t1,dist_y1t3);
                    matrix_cellno_cost(5,:)=([y2,t2,cost_y2t2_y1t1t3]);
                    cost_y2t3_y1t1t2=mycostfunction(area_y2t3,area_y1t1t2,dist_y2t3,dist_y1t1,dist_y1t2);
                    matrix_cellno_cost(6,:)=([y2,t3,cost_y2t3_y1t1t2]);

                    % find best tracking (least cost)
                    matrix_cellno_cost_sort=sortrows(matrix_cellno_cost,3);
                    newyestcellno_nondiv=matrix_cellno_cost_sort(1,1);
                    newtodaycellno_nondiv=matrix_cellno_cost_sort(1,2);
                    
                    NoSwitchHappened=(y2==newyestcellno_nondiv & t3==newtodaycellno_nondiv); %is the old tracking best?
                    if ~NoSwitchHappened
                        allyest=[y1 y2];
                        idxdivP=find(allyest~=newyestcellno_nondiv);
                        newyestcellno_div=allyest(idxdivP);
                        alltoday=[t1 t2 t3];
                        idxdivDE=find(alltoday~=newtodaycellno_nondiv);
                        newtodaycellno_div=alltoday(idxdivDE);
                        fprintf([' \n Found suspicious cell division with large change in area (cutoff ' ...
                            num2str(CUTOFFAREA) ' px): '])
                        fprintf([' \n corrected ' num2str(y1) ' 0 0 ' num2str(t1)])
                        fprintf([' \n           ' num2str(y1) ' 0 0 ' num2str(t2)])
                        fprintf([' \n           ' num2str(y2) ' 0 0 ' num2str(t3)])
                        fprintf([' \n to:       ' num2str(newyestcellno_nondiv) ' 0 0 ' num2str(newtodaycellno_nondiv)])
                        fprintf([' \n           ' num2str(newyestcellno_div) ' 0 0 ' num2str(newtodaycellno_div(1))])
                        fprintf([' \n           ' num2str(newyestcellno_div) ' 0 0 ' num2str(newtodaycellno_div(2))])

                        trackLink(newtodaycellno_nondiv,:)=[newyestcellno_nondiv 0 0 newtodaycellno_nondiv];
                        trackLink(newtodaycellno_div(1),:)=[newyestcellno_div 0 0 newtodaycellno_div(1)];
                        trackLink(newtodaycellno_div(2),:)=[newyestcellno_div 0 0 newtodaycellno_div(2)];
                        
                        % if the previous non-div cell is now dividing:
                        % remove it of the list of non-dividing cells
                        % (don't try new tracking connections)
                        if newyestcellno_div==y2
                            idxnondiv=find(yestcellsStrange_Area_or_Displacement_NonDiv~=newyestcellno_div);
                            yestcellsStrange_Area_or_Displacement_NonDiv=yestcellsStrange_Area_or_Displacement_NonDiv(idxnondiv);
                        end 
                        
                        CellsCorrected=1;
                        break; %allow only for one correction in each cell division -> otherwise the list 
                        % of non-div and div cells would have to be reassigned
                    end
                end % loop over y2=non div cells
                % report , if no correction happened
                if CellsCorrected==0
                    fprintf([' \n Found suspicious cell division with large change in area (cutoff ' ...
                        num2str(CUTOFFAREA) ' px): ' num2str(y1) ...
                        ' -> ' num2str(t1) ' & ' num2str(t2) '. \n Did not correct for it because was ' ...
                        'optimal in cost function'])
                end
                    
            end %loop over y1=strangeDivision
end % if ~isempty(yestCellsstrange.. Div & non-Div)
end % if ~isempty(strange general)
% end % loop SwapCounter

% *****************************************************************************************
% *************** SWAPPED CELLS **************
% *****************************************************************************************

%Brute force:
% Run whole swapping x times so that triplett swaps etc can be detected
% Should be improved by considering a matrix of possible cell swaps, not
% only pairs.

%for Swapcounter=1:p.SwapcounterMax

% -----------------------------------------------------
% Find strange area changes for non-dividing cells
% -----------------------------------------------------
% RECALCULATE STRANGE CELLS (OTHERWISE CORRECTED CELLS IN
% DIVISION_CORRECTION MAY ACCIDENTLY BE IGNORED HERE

% recalculate list because of potential new division associations
%if ~isempty(yestcellsStrangeArea) % has to be done every wround if
%SwapCounter is >1
%    if ~isempty(yestcellsStrangeAreaDiv) % potential need to recalculate which cells are non-div
        % find yesterday cells which are linked to 2 daugthers ->
        % exclude them later
        parent_vec=[];  
        for cellno_maybeparent = all_cellno_yesterday
            % get the cells they link to
            linked_cells = find(trackLink(:,1)==cellno_maybeparent);
            if length(linked_cells) == 2 % two 'today cell' are linked to yesterday cell 'cellno_maybeparent' -> cell division
                parent_vec=[parent_vec; cellno_maybeparent];  
            end
        end
        % get areas of all cells again(today & yesterday)
        allareasyesterday=regionprops(Lc_yesterday_fullsize_resc_centered,'Area'); %use rescaled (?)
        allareasyesterday=[allareasyesterday.Area];
        allareastoday=regionprops(Lc_today_fullsize_centered,'Area');
        allareastoday=[allareastoday.Area];
        %calculate area increase of all tracked cells (if division: accounts for 2 daughters)
        allareas_trackedtotoday=zeros(size(allareasyesterday));  % idx -> cell number of yesterday
        for cellno_yest=all_cellno_yesterday
            linked_cells=find(trackLink(:,1)==cellno_yest);
            allareas_trackedtotoday(cellno_yest)=sum(allareastoday(linked_cells));
            clear linked_cells
        end
        clear cellno_yest
        deltaArea=allareas_trackedtotoday-allareasyesterday;
        %find yesterday cells with strange change in area (e.g.: >150px,<-150px)
        yestcellsStrangeArea=find(abs(deltaArea)>=CUTOFFAREA); %-> yest cell nrs with strange area change

        % find non-dividing(!) yesterday cells which have strange change in area
        % unnecessary here? [NW 2014-05]
        if ~isempty(yestcellsStrangeArea)
            idxStrangeAreaNonDiv=find(~ismember(yestcellsStrangeArea,parent_vec));
            yestcellsStrangeAreaNonDiv=yestcellsStrangeArea(idxStrangeAreaNonDiv);
        end
%    end
%end
% -----------------------------------------------------
% Find Strangely Moving Cells (Displacement (centroid-centroid [not area]))
% -----------------------------------------------------
% 1) if yesterday centroid vector doesn't exist, create it
if ~exist('allcentroidsyesterday')
    allcentroidsyesterday=zeros(max(all_cellno_yesterday),2);
      for cellnorun=all_cellno_yesterday
          cellnoruncentroid=getCentroidsOfCell(Lc_yesterday_fullsize_centered, cellnorun);
          % only 1/2 centroid (and not 1/4 and 3/4) is needed
          cellnoruncentroid=cellnoruncentroid([2 3]); %(y,x)=(row,column)
          allcentroidsyesterday(cellnorun,:)=cellnoruncentroid;
      end
end
% 2) calculate displacements of centroids
allcentroids_trackedtotoday=zeros(size(allcentroidsyesterday));  % idx -> cell number of yesterday
for cellno_yest=all_cellno_yesterday
    linked_cells=find(trackLink(:,1)==cellno_yest);
    % if the cell has dividided take the mean coordinates of all daughter
    % centroids
    if length(linked_cells)==1
        allcentroids_trackedtotoday(cellno_yest,:)=allcentroidstoday(linked_cells,:);
    elseif length(linked_cells)>1
        allcentroids_trackedtotoday(cellno_yest,:)=mean(allcentroidstoday(linked_cells,:));
    end
end
clear cellno_yest
dummy=(allcentroids_trackedtotoday-allcentroidsyesterday);
deltaDisplacement=sqrt(dummy(:,1).^2+dummy(:,2).^2);
%find yesterday cells with strange movement (e.g.: >30px)
% still includes div & non-div cells
yestcellsStrangeDisplacement=find(deltaDisplacement>=CUTOFFDISPLACEMENT); %-> yest cell nrs with strange area change

% % 1) if yesterday centroid vector doesn't exist, create it
%if ~exist('allcentroidsyesterday')
%    allcentroidsyesterday=zeros(max(all_cellno_yesterday),2);
%      for cellnorun=all_cellno_yesterday
%          cellnoruncentroid=getCentroidsOfCell(Lc_yesterday_fullsize_centered, cellnorun);
%          % only 1/2 centroid (and not 1/4 and 3/4) is needed
%          cellnoruncentroid=cellnoruncentroid([2 3]); %(y,x)=(row,column)
%          allcentroidsyesterday(cellnorun,:)=cellnoruncentroid;
%      end
% end


% 2) calculate displacements of centroids
%allcentroids_trackedtotoday=zeros(size(allcentroidsyesterday));  % idx -> cell number of yesterday
%for cellno_yest=all_cellno_yesterday
%    linked_cells=find(trackLink(:,1)==cellno_yest);
%    % if the cell has dividided the centroid thing is pointless -> leave
%    % empty (=0). these dividing cells will be excluded later.
%    if length(linked_cells)==1
%        allcentroids_trackedtotoday(cellno_yest,:)=allcentroidstoday(linked_cells,:);
%    end
%end
%clear cellno_yest
%dummy=(allcentroids_trackedtotoday-allcentroidsyesterday);
%deltaDisplacement=sqrt(dummy(:,1).^2+dummy(:,2).^2);
% %find yesterday cells with strange movement (e.g.: >30px)
% %  this also includes thedividing cells because the daughter centroid is set
% % to (0,0). they will be excluded later.
%yestcellsStrangeDisplacement=find(deltaDisplacement>=cutoff_Displacement); %-> yest cell nrs with strange area change



% --------------------------------------------------
% ------- combine strange yest. cells and restrict to non-dividing cells ---
% --------------------------------------------------
yestcellsStrange_Area_or_Displacement=unique([yestcellsStrangeArea, yestcellsStrangeDisplacement']);

if ~isempty(yestcellsStrange_Area_or_Displacement)
    idxStrangeNonDiv=find(~ismember(yestcellsStrange_Area_or_Displacement,parent_vec));
    yestcellsStrange_Area_or_Displacement_NonDiv=yestcellsStrange_Area_or_Displacement(idxStrangeNonDiv);
end



% -----------------------------------------------------
%  ----- Attempt correction of swapped cells --------
% -----------------------------------------------------

if ~isempty(yestcellsStrange_Area_or_Displacement)
    % Non-Div cells move strange or change area a lot:
    if ~isempty(yestcellsStrange_Area_or_Displacement_NonDiv)

        % 6) loop over all non-dividing strange area cells -> check for swaps
        % create cost function for possible tracking links. THIS IS VERY
            % EMPIRICAL AND MIGHT HAVE TO BE ADJUSTED!!
            % also no renormalization of areas & displacements
            % Normalization Parameters: defined above
% MAYBE TODO: EXTRA EMPHASIS ON THAT THE FIRST TRACK (LOWER CELL NO) IS CORRECT
% (the higher nr can be corrected again in a subsequent loop)
            % **** IMPORTANT ****
            % First yest cell (y1, lower cenn nr counter): deltaarea1,
            % deltadist1
            % KEEP CORRECT ORDER OF INPUT ARGUMENTS!
            % First cell is weighed stronger (2nd cell has 2nd chance to
            % retrack
            % ******************
            mycostfunction = @ (deltaarea1, deltaarea2, deltadist1, deltadist2) ...
               (deltaarea1 + deltaarea2)/NORMALAREACHANGENONDIV + (deltadist1+deltadist2)/NORMALDISPLACEMENTNONDIV ...
                + deltaarea1/NORMALAREACHANGENONDIV *3 + deltadist1/NORMALDISPLACEMENTNONDIV *3; %overemphasize y1 cell (*3)    
           %+ deltadist1/normalDisplacementNonDiv *10; %overemphasize dist y1 cell (*10)
               
           
            
           %blubbmycostfunction = @ (deltaarea1, deltaarea2, deltadist1, deltadist2) ...
           %    (deltaarea1 + deltaarea2)/normalAreaChangeNonDiv + (deltadist1+deltadist2)/normalDisplacementNonDiv;
 
            % input: all area changes (2), all displacements (2)
            % output: cost of tracking connection

            % cells which were already relinked
            yest_cells_to_skip=[];

            for idx_y1cell=1:length(yestcellsStrange_Area_or_Displacement_NonDiv)
                y1=yestcellsStrange_Area_or_Displacement_NonDiv(idx_y1cell);
                if ismember(y1, yest_cells_to_skip)
                    continue
                end
                CellsCorrected=0;
                t1=find(trackLink(:,1)==y1);
                    % loop over all strange-area non-dividing cells (avoid
                    % double checking)
                    for idx_y2cell=idx_y1cell+1:length(yestcellsStrange_Area_or_Displacement_NonDiv)
                        y2=yestcellsStrange_Area_or_Displacement_NonDiv(idx_y2cell);
                        %disp('...')
                        %disp(['y1 cell: ' num2str(y1) ' y2 cell: ' num2str(y2)])
                        % get all cells of today which are linked to y1 & y2
                        t2=find(trackLink(:,1)==y2);
                        if length(t1)~=1 | length(t2)~=1
                            fprintf(' \n Warning: Error in tracking (area) correction of non-dividing cells.')
                            fprintf('\n y1 & y2 should have one link. Will continue with next y1&y2. Barren cell?')
                            continue

                        end
                        % permute over tracking possibilities and find all abs. area
                        % changes & centroid displacements (could be improved with
                        % a loop/matrix)
                        % for convenience and readability, precalculate all
                        % possible abs area changes and centroid<->area
                        % displacements
                        area_y1t1=abs(allareastoday(t1)-allareasyesterday(y1));   % area y1->t1 (no division)
                        area_y1t2=abs(allareastoday(t2)-allareasyesterday(y1));
                        area_y2t1=abs(allareastoday(t1)-allareasyesterday(y2));
                        area_y2t2=abs(allareastoday(t2)-allareasyesterday(y2));

                        dist_y1t1=GetDistance_Point_from_Area(allpixellistYesterday(y1).PixelList(:,2), allpixellistYesterday(y1).PixelList(:,1), allcentroidstoday(t1,1), allcentroidstoday(t1,2));
                        dist_y1t2=GetDistance_Point_from_Area(allpixellistYesterday(y1).PixelList(:,2), allpixellistYesterday(y1).PixelList(:,1), allcentroidstoday(t2,1), allcentroidstoday(t2,2));
                        dist_y2t1=GetDistance_Point_from_Area(allpixellistYesterday(y2).PixelList(:,2), allpixellistYesterday(y2).PixelList(:,1), allcentroidstoday(t1,1), allcentroidstoday(t1,2));
                        dist_y2t2=GetDistance_Point_from_Area(allpixellistYesterday(y2).PixelList(:,2), allpixellistYesterday(y2).PixelList(:,1), allcentroidstoday(t2,1), allcentroidstoday(t2,2));

                        % calculate cost functions of tracking & create an array to easily
                        % recover the cell numbers
                        matrix_cellno_cost=zeros(2,3); % per row: first non-div yest cell  --- non-div today cell --- cost
                        % 2 options: y1->t1 or y1->t2
                        
                        % *********** always begin with y1 area & y1 dist !!! ***********
                        cost_y1t1_y2t2=mycostfunction(area_y1t1,area_y2t2,dist_y1t1,dist_y2t2); % y1->t1, y2->t2
                        matrix_cellno_cost(1,:)=([y1,t1,cost_y1t1_y2t2]);
                        cost_y1t2_y2t1=mycostfunction(area_y1t2,area_y2t1,dist_y1t2,dist_y2t1);
                        matrix_cellno_cost(2,:)=([y1,t2,cost_y1t2_y2t1]);

                        % find best tracking (least cost)
                        matrix_cellno_cost_sort=sortrows(matrix_cellno_cost,3);
                        newyestcellno_1=matrix_cellno_cost_sort(1,1);
                        newtodaycellno_1=matrix_cellno_cost_sort(1,2);   
                        NoSwitchHappened=(y1==newyestcellno_1 & t1==newtodaycellno_1) | (y2==newyestcellno_1 & t2==newtodaycellno_1);
                        %disp(matrix_cellno_cost)
                        if ~NoSwitchHappened
                            allyest=[y1 y2];
                            idxy=find(allyest~=newyestcellno_1);
                            newyestcellno_2=allyest(idxy);
                            alltoday=[t1 t2];
                            idxt=find(alltoday~=newtodaycellno_1);
                            newtodaycellno_2=alltoday(idxt);
                            fprintf([' \n Found suspicious change in cell area (>' ...
                                num2str(CUTOFFAREA) ' px) or displacement (>' num2str(CUTOFFDISPLACEMENT) ...
                                 ' px) : '])
                            fprintf([' \n corrected ' num2str(y1) ' 0 0 ' num2str(t1)])
                            fprintf([' \n           ' num2str(y2) ' 0 0 ' num2str(t2)])
                            fprintf([' \n to:       ' num2str(newyestcellno_1) ' 0 0 ' num2str(newtodaycellno_1)])
                            fprintf([' \n           ' num2str(newyestcellno_2) ' 0 0 ' num2str(newtodaycellno_2)])

                            CellsCorrected=1;

                            trackLink(newtodaycellno_1,:)=[newyestcellno_1 0 0 newtodaycellno_1];
                            trackLink(newtodaycellno_2,:)=[newyestcellno_2 0 0 newtodaycellno_2];

                            % if switch happened move on to next y1 cell
                            % (yesterday). only allow one switch to happen. That is with
                            % the first y2 cell where the swap is better than
                            % the original trackings. Note: This could be
                            % improved by finding global best switch, but
                            % typically only one switch should be better than
                            % the current tracking
                            % remove the swapped y2 cell from the investigated
                            % list:
                            yest_cells_to_skip=[yest_cells_to_skip, y2]; % necessary at all? NW 2014-05
                            break

                        end
                    end % loop over y2=non div cells
                    % if cell was not corrected, report this:
                    if CellsCorrected==0 & ~isempty(idx_y2cell)
                        fprintf([' \n Found suspicious change in cell area (>' ...
                                num2str(CUTOFFAREA) ' px) or displacement (>' num2str(CUTOFFDISPLACEMENT) ...
                                ' px): ' num2str(y1) ...
                                ' -> ' num2str(t1) '. \n Did not correct for it because was ' ...
                                'optimal in cost function'])
                    end                      
                end %loop over y1=non div cells
    end % loop over yestcellsStrange..nonDiv
    
end % if ~isemptys(strange cells..)


end % loop SwapCounter
% -------------------------------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % CORRECT FORMAT OF DIVIDING CELLS IN TRACKING
  %------------------------------------------------------------------------
  % now loop over all cellno in yesterday, if any occurs twice, divide
  for cellno = all_cellno_yesterday
    if cellno==0, continue; end % background
    
    % how often does this cellno come back?
    cellno_children = find(trackLink(:,1)==cellno);
    nr_occur = length(cellno_children);
    
    % if cellno occurs 2 times, cell has divided, and trackLink should be adjusted
    if nr_occur==2
      % need to divide
      cellno_today_1 = trackLink(cellno_children(1),4);
      cellno_today_2 = trackLink(cellno_children(2),4);
      trackLink(cellno_children(1),:) = [0 cellno 0 cellno_today_1];
      trackLink(cellno_children(2),:) = [0 0 cellno cellno_today_2];
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
% Convert matching results to schnitzcells-format lineage
MW_makeSchnitzFileFromTracking(p, opts);
%--------------------------------------------------------------------------


% **************************************************************************************************
% **********************************END MAIN FUNCTION **********************************************
% **************************************************************************************************

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
% write Lc into full size image
  % first, check if coordinates are out of bounds -> don't rescale image
  minrow=min(fy)+offset_y;
  maxrow=max(fy)+offset_y;
  mincol=min(fx)+offset_x;
  maxcol=max(fx)+offset_x;
  sizeLcfull=size(Lc_fullsize);
  if minrow<1 | maxrow>sizeLcfull(1) | mincol<1 | maxcol>sizeLcfull(2)
      fprintf('Image too large for centering (idx out of bounds). Will use non-centered image.')
      Lc_fullsize_centered=Lc_fullsize;
  else
      Lc_fullsize_centered( minrow:maxrow, mincol:maxcol ) = Lc_fullsize( min(fy):max(fy), min(fx):max(fx) );
  end
    
  DEBUGLOADLCDATA=0;
  if DEBUGLOADLCDATA
     figure
     imagesc(Lc_fullsize_centered)
     title('full size centered')
     figure
     imagesc(Lc_fullsize)
     title('full size. no centering')
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rescale the yesterday-image to have the same size as the today image (to
% acocunt for cell growth).
% recenter it.
function [YestIm_fullsize_resc]=rescaleLcData(YestIm_orig,TodayIm_orig,YestIm_fullsize)
  % determine rescale factor (mass increase)
  totalMass_yesterday=sum(sum(YestIm_orig>0));
  totalMass_today=sum(sum(TodayIm_orig>0));
  rescaleLc=totalMass_today/totalMass_yesterday;
  YestIm_resc=imresize(YestIm_orig,rescaleLc,'nearest');
  % initialize full size image
  sizeLcfull=size(YestIm_fullsize);
  YestIm_fullsize_resc=zeros(sizeLcfull);
  % get weighted center of cells in Lc
  [fy, fx]= find(YestIm_resc>0);
  center_cells_x = round( mean(fx) ); 
  center_cells_y = round( mean(fy) ); 

  % determine offset
  offset_x = round( size(YestIm_fullsize_resc,2)/2 - center_cells_x);
  offset_y = round( size(YestIm_fullsize_resc,1)/2 - center_cells_y);
  
  % write Lc into full size image
  % first, check if coordinates are out of bounds -> don't rescale image
  minrow=min(fy)+offset_y;
  maxrow=max(fy)+offset_y;
  mincol=min(fx)+offset_x;
  maxcol=max(fx)+offset_x;
  if minrow<1 | maxrow>sizeLcfull(1) | mincol<1 | maxcol>sizeLcfull(2)
      fprintf('Image too large (idx out of bounds for rescaling). Will use non-rescaled image for yesterday-image.')
      YestIm_fullsize_resc=YestIm_fullsize;
  else
      YestIm_fullsize_resc( minrow:maxrow,mincol:maxcol) = YestIm_resc( min(fy):max(fy), min(fx):max(fx) );
  end
    
  DEBUGRESCALE=0;
  if DEBUGRESCALE
     figure
     imagesc(YestIm_fullsize_resc)
     title(['yesterday fullsize. rescaled by factor ' num2str(rescaleLc)])
     figure
     imagesc(YestIm_fullsize)
     title('yesterday fullsize. no rescale')
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
xmin= max(min(fx) - extra, 1); %rows (so in std understanding the y-axis, comment NW2013-12)
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
  centroids(2:3) = [(px(cen(2)) + xmin - 1) (py(cen(2)) + ymin - 1)]; %NW2013-12:careful:
  centroids(4:5) = [(px(cen(1)) + xmin - 1) (py(cen(1)) + ymin - 1)]; % here: x=row, y=column
  centroids(6:7) = [(px(cen(3)) + xmin - 1) (py(cen(3)) + ymin - 1)];
else
  disp('Thin to short?');
  centroids = [cellno round(mean(fx)) round(mean(fy)) round(mean(fx)) round(mean(fy)) round(mean(fx)) round(mean(fy))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the minimal distance of one point (typically the backbone
% centroid) from an area
 function mycelldist=GetDistance_Point_from_Area(myarearow, myareacol, mycentroidrow, mycentroidcol)
            distvec=sqrt((myarearow-mycentroidrow).^2+(myareacol-mycentroidcol).^2);
            mycelldist=min(distvec); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

