function []=NW_makeMovie_branchGroups(p,branchGroups,varargin)
% ***************************************************
% follows strongly Daan's makeMovie 
% removal of many options. cells will automatically be colored by number of
% the branchGroup they are in.
% Legend displays which color belongs to which branchgroup and
% start-schnitz
% ***************************************************
%
% NW_makeMovie_branchGroups makes movie of cells. Gets data from segmentation files and 
% schnitzcells and saves images colored accoding to branchGroup.
% Full tracking is not required, but 'singleCell' tracking is.
% 
%
% Works only for up to 12 branchGroups so far! (otherwise change colorMap)
%
% It might be useful to label cells with branchnumber and not schnitznumber
%
% If cells are not part of branch groups (i.e. out of fitting time or
% removed), they will appear grey
%
% It is important that the colormap in
% DJK_plotcrosscorrelation_standard_error_store is the same!
%
%
% REQUIRED ARGUMENTS:
%
% p
% branchGroups        output of DJK_divide_branch_data
%
% OPTIONAL ARGUMENTS:
%
% 'manualRange'       Allows to analyze a subset of frames
%
% 'onScreen' = 1      Shows images onScreen, instead of writing to disk.
%
% 'stabilize' = 0     Turns of stabilize image (less shaking of cells) default:true
%
% 'DJK_saveDir'       Directory where images should be saved. Defaults to
%                     "p.analysisDir \ movies \ movieName \"
%
% 'schnitzcells'      Use the provided as schnitzcells or lineage data
%
% 'problemCells'      array with problematic cells (+ frame nr), so that
%                     only schnitz nrs will be written in these.
%                     Still in program, but is not supposed to be used in
%                     the branchGroup version (supposedly works
%                     nevertheless)


% ****** pre-set options (to stay close to Daan's version) **********
mode='branchgroupnr';   % Color of cell is determined by number of branchGroup
wrNum='schnitz';        % Writes schnitz number in cell when cell first occurs
case_writeNum=1;        % add numbers to cells
case_writeSchnitz = true; % add Schnitz numbers
case_writeSchnitzAll= true;  % no idea for what...(NW)

% ******************************************************************





%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2;
functionName = 'NW_makeMovie_branchGroups';

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

% Get names of segmentation files in segmentation directory
outprefix = [p.movieName 'seg'];
D = dir([p.segmentationDir, outprefix, '*.mat']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;

% If problemCells mode, numbers are only written for cells/frames in provided problemCells
if ~existfield(p,'problemCells')
  p.problemCells = false;
end

  
% If explicit manualRange is not given, take all segmentation files
if ~existfield(p,'manualRange')
  segNameStrings = char(S);
  p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
  if p.problemCells % problemcells, select frames with problem cells, also before & after
    errorFramesRange = [];
    errorFrames = unique(p.problemCells(:,2))';
    for fr = errorFrames
      idx = find(p.manualRange == fr);
      if (idx-1>=1)                     , errorFramesRange = [errorFramesRange p.manualRange(idx-1)]; end
      if (idx)                          , errorFramesRange = [errorFramesRange p.manualRange(idx)]; end
      if (idx+1<=length(p.manualRange)) , errorFramesRange = [errorFramesRange p.manualRange(idx+1)]; end
    end
    p.manualRange = unique( errorFramesRange );
%     fr_min = p.manualRange(1); fr_max = p.manualRange(end); p.manualRange = [];
%     errorFrames = unique(p.problemCells(:,2))';
%     for fr = errorFrames
%       p.manualRange = [p.manualRange fr-1 fr fr+1];
%     end
%     p.manualRange = unique(  p.manualRange(p.manualRange>fr_min-1 & p.manualRange<fr_max+1)  );
  end
end
disp(['Analyzing ' num2str(length(p.manualRange)) ' frames ', num2str(p.manualRange(1)), ' to ', num2str(p.manualRange(end))]);

% Check whether there are frames without a corrected segmentation
for frameNum = p.manualRange
  clear Lc;
  load([p.segmentationDir, p.movieName, 'seg', str3(frameNum)]);
  if ~(exist('Lc') == 1) 
    disp(['Not corrected segmentation in seg file : ', p.movieName, 'seg', str3(frameNum)]);
  end
end


% If onScreen, nothing is saved to disc
if ~existfield(p,'onScreen')
  p.onScreen = false;
end

% If stabilize, will try to stabilize image
if ~existfield(p,'stabilize')
  p.stabilize = true;
end


% Set lineage file to load.
% ********** _lin.mat is not enough here since you need time info (->
% fluorinfo) ********** (NW 2012-03)
if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end

% If schnitzcells / lineage file is provided, copy to lineage and remove from p
if existfield(p,'schnitzcells')
  lineage.schnitzcells = p.schnitzcells;
  p = rmfield(p, 'schnitzcells');
end

% saveDirectory and Name
filenameBase=['BranchGroupMovies_SchnitzNum'];
if (p.stabilize)
  filenameBase = [filenameBase '_Stabile'];
end
if (p.problemCells)
  filenameBase = [filenameBase 'Problem'];
end

% If explicit DJK_saveDir is not given, define it
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'branches' filesep  filenameBase filesep];
end

% Make sure DJK_saveDir exists, else make it
if isfield(p,'DJK_saveDir') & ~p.onScreen
  % make sure every directory field has a trailing filesep
  pathlength = length(p.DJK_saveDir);
  if (p.DJK_saveDir(pathlength) ~= filesep)
    p.DJK_saveDir = [p.DJK_saveDir filesep];
  end
  % if directory doesn't exist, create it
  if exist(p.DJK_saveDir)~=7
    [status,msg,id] = mymkdir([p.DJK_saveDir]);
    if status == 0
      disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
      return;
    end
  end
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PREPARING TRACKING / SCHNITZCELLS FILE, SLOOKUP AND COLORMAP
%--------------------------------------------------------------------------
if ~exist('lineage') % if schnitzcells provided, no need to load
    % load lineage
    if ~(exist(p.schnitzName)==2)
        error(['Could not read lineage file ''' p.schnitzName '''.']);
    end
    disp(['Loading lineage from file ''', p.schnitzName, '''.']);
    lineage = load(p.schnitzName);  
end


% slookup gives the schnitznr for each cellno in each frame.
slookup = makeslookup(lineage.schnitzcells);


% findstartschnitz associates the starting/parent schnitz of the
% branchGroup to each schnitznumber
% if schnitz is too old to appear in BranchGroup the starting schnitz will
% be associated to =0
startschnitz_array=findstartschnitz(lineage.schnitzcells,branchGroups);

% take same colorMap as in DJK_plot_crosscorrelation_standarderror_store
% (NW 2012/03))
    myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];%myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
% shift by one to let schnitzes that are not in branchGroup appear grey ('0' in startschnitz_array)
myColor=[0.5 0.5 0.5; myColor];

% convert to Daan's version
mymap=[0 0 0;myColor];% first entry for background color
number_of_colors=size(mymap,1);
% explicitly set some special colors in the colormap
color_background  = 1;
color_text        = number_of_colors+1;
legend_background = number_of_colors+2;
color_legendtext        = number_of_colors+3;

mymap(color_background,:) = [0   0   0  ];
mymap(color_text,:)       = [1   1   1  ];
mymap(legend_background,:)       = [1   1   1  ];
mymap(color_legendtext,:)=[0 0 0];

%--------------------------------------------------------------------------
% creates n x 2 array with frame number and associated time
%--------------------------------------------------------------------------
% (allows to check later if imaged frame is part of branchgroup range)

% Takes time points of first existent fluor-color (usually fluor1).

%get existent fluor color
if strcmp(p.fluor1,'none')==0
    fluor_time=[upper(p.fluor1) '_time'];
    fluor_frame=[upper(p.fluor1) '_frames'];
    disp(['Using ' p.fluor1 ' frames to get time']);
elseif strcmp(p.fluor2,'none')==0
    fluor_time=[upper(p.fluor2) '_time'];
    fluor_frame=[upper(p.fluor2) '_frames'];
    disp(['Using ' p.fluor2 ' frames to get time']);
elseif strcmp(p.fluor3,'none')==0
    fluor_time=[upper(p.fluor3) '_time'];
    fluor_frame=[upper(p.fluor3) '_frames'];
    disp(['Using ' p.fluor3 ' frames to get time']);
else
    disp(['You don''t have any fluor colors. Don''t know to which time field to use. Exiting...']);
    return
end
% matrix with [frames, real time]
timetable=associate_frames_realtime(lineage.schnitzcells,fluor_time, fluor_frame);
%get time range of branchgroup-data
branchtime_min=branchGroups(1).branches(1).(fluor_time)(1);
branchtime_max=branchGroups(1).branches(1).(fluor_time)(length(branchGroups(1).branches(1).(fluor_time)));
%add column to timetable which determines if frame is within
%branchgroup time range
inrange=(timetable(:,2)>=branchtime_min & timetable(:,2)<=branchtime_max);
%shift by 1: wrong frame association in schnitzcells (comment
%NW2012-03): when schnitz has 420 as first entry in frames, it
 %actually already appears in frame 419!!
timetable(:,1)=timetable(:,1)-1;
timetable=[timetable,inrange]
%find min and max frame of branchgroup
for i=1:size(timetable,1)
    if timetable(i,3)==1;
        framemin=timetable(i,1)
        break
    end
end
for i=size(timetable,1):(-1):1
    if timetable(i,3)==1;
        framemax=timetable(i,1)
        break
    end
end
%--------------------------------------------------------------------------
% LOOPING OVER SEGMENTATION FILES
%--------------------------------------------------------------------------
% loop over the frames
for fr = p.manualRange
  % load the seg file
  clear Lc LNsub creg yreg savelist rect newrect oldrect;
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load(name);

  % if segmentattion is not approved, get unapproved segmentation 
  if ~(exist('Lc') == 1) 
    Lc = LNsub;
  end

  % get the unique nrs in segmentation (u)
  u = unique(Lc); 
  % get 0 out of the list of unique cellno
  u = setdiff(u,0);

  % get the segmented image (subim) 
  subim = double(Lc);
  % subim_color will be the colored image
  subim_color = zeros(size(subim)); % 0 will be background color
  % If number has to be written keep a separate image to put them into
  if (case_writeNum)
    subim_text = zeros(size(subim));
  end


  %------------------------------------------------------------------------
  % LOOPING OVER CELL NRS
  %------------------------------------------------------------------------
  for uind = 1:length(u),
    % cellno uind in framenum fr is 'just born', 'growing' or 'dividing'
    temp_schnitz = slookup(fr+1,u(uind));
    if temp_schnitz==0
      disp(['In frame ' str3(fr) ' cellno ' str3(uind) ' is not linked to schnitz']);
      continue;
    end
    temp_schnitz_P = lineage.schnitzcells(temp_schnitz).P;
    temp_schnitz_frames = lineage.schnitzcells(temp_schnitz).frames;
    temp_schnitz_cellno = lineage.schnitzcells(temp_schnitz).cellno;
    frIndex = find(temp_schnitz_frames==(fr+1));

    % fl are the locations in image where cellno uind is located
    fl = find(subim==u(uind));

    %----------------------------------------------------------------------
    % DO THE COLORING
    %----------------------------------------------------------------------
    
    %get Schnitznumber
    myschnitz=slookup(fr+1,u(uind));
    %convert to branchgroupindex
    idx=find(myschnitz==startschnitz_array(:,2));
    % color according to branchGroupNr+1 (+2 cause else schnitzes without branchGroup-parent are black (and index starts at =0))
    subim_color(fl) = startschnitz_array(idx)+2;   
   
    %----------------------------------------------------------------------


    
    %----------------------------------------------------------------------
    % WRITE NUMBERS (SCHNITZ NUMBERS OR CELL NUMBERS)
    %----------------------------------------------------------------------
    if (case_writeNum)
      ref_fr = find(temp_schnitz_frames==fr+1);

      % if we are in problem mode check whether correct frame/cell, else continue 
      if p.problemCells
        errorFramesForThisSchnitz = p.problemCells(find(p.problemCells(:,1)==temp_schnitz),2);
        if length(intersect( fr, errorFramesForThisSchnitz)) | ~case_writeSchnitz % if cellno, write all of them
        else
          continue;
        end
      end
      
      if (case_writeSchnitz) % write schnitz numbers
        nr_im = DJK_getImageOfNr( temp_schnitz );
      else % write cellno in each cell
        nr_im = DJK_getImageOfNr( temp_schnitz_cellno(ref_fr) );
      end
      
      % coordinates
      nr_x = lineage.schnitzcells(temp_schnitz).cenx(ref_fr);
      nr_y = lineage.schnitzcells(temp_schnitz).ceny(ref_fr);

      % write number
      if (temp_schnitz_frames(1)==fr+1 | case_writeSchnitzAll | ~case_writeSchnitz)
        subim_text = DJK_imagePlace(subim_text, nr_im, nr_x, nr_y);
      end
    end
    %----------------------------------------------------------------------
    
  end
  %--------------------- END LOOPING OVER CELL NRS  -----------------------

  
  %------------------------------------------------------------------------
  % FINISH IMAGE
  %------------------------------------------------------------------------
    if (case_writeNum) % put schnitz/cellno numbers in image
      subim_color(find(subim_text>0)) = color_text;
    end
    
    if (p.onScreen)
      % create legend
      mylegend=['''branchgroup 1. startschnitz ' num2str(branchGroups(1).parent_cell) ''''];
      for i=2:length(branchGroups)
          mylegend=[mylegend, ', ''branchgroup ', num2str(i) '. startschnitz ' , num2str(branchGroups(i).parent_cell),''''];
      end
              
      figure;
      image(ind2rgb(subim_color,mymap));
      axis off;
      axis image;
      %dumb solution to get sth to write a legend on
      hold on
      for i=1:length(branchGroups)
          plot(1,1,'Color',mymap(i+2,:))
      end
      eval(['legend(' mylegend ',''Location'',''SouthEastOutside'')']);
    else
      % enlarge to original size
      im_color = DJK_enlargeImage(subim_color, phaseFullSize, rect, p.stabilize);
      
      % manually add legend (burn into image)
      % white background
      backXmax=size(im_color,2)-20; backYmax=size(im_color,1)-20;sizeX=300; sizeY=30*length(branchGroups);
      backXcent=backXmax-0.5*sizeX;backYmin=backYmax-sizeY; backXmin=backXmax-sizeX;% for convenience
      % *** cheating start ***     
      %backYmin=10; %crude forcing of legend coordinnates to not decrease
      %sizeY=backYmax-10; %below 0 - if many branches exist
      % *** cheating stop ***
      im_color(backYmax-sizeY:backYmax,backXmax-sizeX:backXmax)=legend_background;
      % add legend text
      for i=1:length(branchGroups)
          clear legendImText
          legendImText=text2im(['branchgr ', num2str(i) '. schnitz ' , num2str(branchGroups(i).parent_cell)]); %get image of text
          legendImText=imresize(legendImText,0.7,'nearest');
          legendImText=legendImText*legend_background+(1-legendImText)*color_legendtext; % adjust values to get right color
          im_color=DJK_imagePlace(im_color,legendImText,backXcent+10,backYmin+10+(i-1)*30);
      end
      %draw colored squares according to cell colors
      for i=1:length(branchGroups)
        im_color(backYmin+5+(i-1)*30:backYmin+15+(i-1)*30,backXmin+5:backXmin+15)=i+2;
      end
      
      % check whether frame number is used in fitRange of BranchGroups. 
      % Otherwise add warning text
      if fr<framemin | fr>framemax
          %shift by 1: wrong frame association in schnitzcells (comment
          %NW2012-03): when schnitz has 420 as first entry in frames, it
          %actually already appears in frame 419!!
          outOfRangeImText=text2im('frame not in branchgroup time');
          outOfRangeImText=imresize(outOfRangeImText,0.7,'nearest');
          outOfRangeImText=outOfRangeImText*color_background+(1-outOfRangeImText)*color_text; % adjust values to get right color
          im_color=DJK_imagePlace(im_color,outOfRangeImText,250,50);
      end
      
      % save image
      imwrite(ind2rgb(im_color,mymap),[p.DJK_saveDir filenameBase str3(fr) '.png'],'png');
     
      disp(['Writing image : ', [filenameBase str3(fr) '.png']]);
    end

 
end
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makesloopup makes array: schnitz = slookup( frame+1, cellno)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slookup = makeslookup(s);
for i = 1:length(s),
  for j = 1:length(s(i).frames),
    if ~isempty(s(i).frames),
      if s(i).frames
        slookup(s(i).frames(j),s(i).cellno(j))=i;
      end;
    end;
  end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% branchGroup to each schnitznumber
% if schnitz is too old to appear in BranchGroup the starting schnitz will
% be associated to =0
% e.g. [0  1; 0  2; 0  3; 1  4; ....; 1  8; 1  9; 2  10; ... ;
% branchgroupnr schnitznr; ...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startschnitz_array=findstartschnitz(myschnitzcells,branchGroups);
startschnitz_array=zeros(length(myschnitzcells),2);
for i=1:length(myschnitzcells)
   startschnitz_array(i,2)=i;
   for grouprun=1:length(branchGroups)
       for branchrun=1:length(branchGroups(grouprun).branches)
           clear idx;
           idx=find(i==branchGroups(grouprun).branches(branchrun).schnitzNrs);
           if ~isempty(idx) % schnitz appears in this branch
               startschnitz_array(i,1)=grouprun; % get number of branchgroup
               continue
           end
       end
       if ~isempty(idx), continue, end; % jump to next schnitz
   end
   if ~isempty(idx), continue, end; % jump to next schnitz
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time point to each fluor frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timetable=associate_frames_realtime(myschnitzcells,myfluortime,myfluorframe);
myframes=[];
mytime=[];
for i=1:length(myschnitzcells)
    myframes=[myframes; myschnitzcells(i).(myfluorframe)'];
    mytime=[mytime; myschnitzcells(i).(myfluortime)'];
end
myframes=unique(myframes);
mytime=unique(mytime);
timetable=[myframes, mytime];

        


