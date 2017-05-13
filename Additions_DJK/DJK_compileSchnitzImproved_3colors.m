% DJK_compileSchnitzImproved_3colors extracts data out of lineage & segmentation files.
%
% NW2012-08: adapted such that program can also cope with different
% times/frequencies for different fluorescence images (cf focus test)
%
% Copied / Rewritten from COMPILESCHNITZ.
%
% It adds up to three fluorescence colors (depending on whether colors are
% specified in 'p')
%
% OUTPUT
% 'p'   
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'p'
%
% OPTIONAL ARGUMENTS:
% 'schnitzName'       name of annotated output schnitz file (default: [movieName '-Schnitz.mat'])
% 'lineageName'       name of un-annotated input lineage name (default: [movieName '_lin.mat'])
% 'quickMode'=1       will run in quick mode, i.e. load the previously made schnitzcells file. (default = 0)
% 'micronsPerPixel'   default = 0.04065 (1.5x magnification of Killer Mike)
% 'cropLeftTop'

function [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,varargin) 

%%
%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_compileSchnitzImproved_3colors';

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

%%
%--------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% Set default parameter values if they don't exist yet
if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end
if ~existfield(p,'quickMode')
  p.quickMode = 0;
end
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end
if ~existfield(p,'micronsPerPixel')
  p.micronsPerPixel = 0.04065; % CoolSnap camera, setup 1
  disp('WARNING: p.micronsPerPixel was not set. Assuming CoolSnap camera, setup 1.');
end
if ~existfield(p,'cropLeftTop')
  % if not set
  disp('Warning! p.cropLeftTop not set');
  cropLeftTop = [1,1];
elseif isempty(p.cropLeftTop)
  % if empty (which also means no cropping was applied)
  cropLeftTop = [1,1];
else
  % if value was given
  cropLeftTop = p.cropLeftTop;
end
%--------------------------------------------------------------------------
  

%--------------------------------------------------------------------------
% QuickMode, so only load previously extracted data
%--------------------------------------------------------------------------
if p.quickMode 
  disp('quick mode on! skipping data extraction! Loading data from file');
  % Loading an existing schnitz file which contains the image-derived fields
  if exist(p.schnitzName)~=2
    error(['Could not read ' p.schnitzName ' , which is required for quick mode']);
  else
    load(p.schnitzName);
    disp(['Load from ''' p.schnitzName ''' completed...']);
  end
  return; % exit here
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOADING OF SCHNITZCELLS LINEAGE AND CREATION OF lincellnum
%--------------------------------------------------------------------------
% Will use a lincellnum structure, where for only used frames (not necessarly
% starting with image 001 (frame 2), the link between schnitznum and cellno
% is saved: lincellnum {a} (b) returns schnitznum of cellno b in frame a

% check for and load lineage file. This will contain 'schnitzcell' with basic data in it 
if ~(exist(p.lineageName)==2)
  error(['Could not read lineage file ''' p.lineageName '''.']);
end
disp(['Loading lineage from file ''' p.lineageName ''' and extracting image-based fields...']);
load(p.lineageName);

% Remove some old obsolete fields
schnitzcells = rmfield(schnitzcells, ['ang      '; 'len      ']); % DJK 090619 'cenx     '; 'ceny     '; 

% Get trackRange (frames that will be extracted). Note: image 001 is frame 2 -> -1
trackRange = sort(unique([schnitzcells.frame_nrs])); % MW 2014/06/11 removal N+1 bug

% initialize lincellnum to have zero'd arrays for each frame 
lincellnum = {};
for lincellnumIndex = 1:length(trackRange)
  segdata = load([p.segmentationDir,p.movieName,'seg',str3(trackRange(lincellnumIndex))]);
  if ~isfield(segdata,'Lc');
    segdata.Lc = segdata.LNsub;
  end
  lincellnum{lincellnumIndex} = zeros([1 max2(segdata.Lc)]);
end

% Now, step through each schnitz, store schnitz number in lincellnum
for schnitznum = 1:length(schnitzcells)
  s = schnitzcells(schnitznum);
  for age = 1:length(s.frame_nrs) % MW 2014/06/11 removal N+1 bug
    framenum = s.frame_nrs(age); % MW 2014/06/11 removal N+1 bug

    % look up that frame number's index within trackRange
    lincellnumIndex = find(trackRange==framenum);
    cellnum  = s.cellno(age);
    lincellnum{lincellnumIndex}(cellnum) = schnitznum;
  end
end
%--------------------------------------------------------------------------


% -------------------------------------------------
% create variable names
% -------------------------------------------------
%%
reg1=genvarname([p.fluor1 'reg']);  % creates 'nonereg' if fluorcolor not existent
reg2=genvarname([p.fluor2 'reg']);  % -> distinguish later if real fluorcolor
reg3=genvarname([p.fluor3 'reg']);  %
fluor1_frames_all=genvarname([upper(p.fluor1) '_frames_all']);  %'NONE_frames_all' if color non-existent
fluor2_frames_all=genvarname([upper(p.fluor2) '_frames_all']);  % e.g. 'Y_frames_all' if yellow fluor
fluor3_frames_all=genvarname([upper(p.fluor3) '_frames_all']);
fluor1_frames=genvarname([upper(p.fluor1) '_frames']);
fluor2_frames=genvarname([upper(p.fluor2) '_frames']);
fluor3_frames=genvarname([upper(p.fluor3) '_frames']);
fluor1_time=genvarname([upper(p.fluor1) '_time']);
fluor2_time=genvarname([upper(p.fluor2) '_time']);
fluor3_time=genvarname([upper(p.fluor3) '_time']);

% --------------------------------------------------

% --------------------------------------------------
% create control variables for check if right fluorescence colors were
% added to 'p'
% --------------------------------------------------
fluor1counter=0; fluor2counter=0; fluor3counter=0;
% --------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER FRAMES IN TRACKRANGE
%--------------------------------------------------------------------------
%%
disp(['-----------------------------------------------------------']);
for i = 1:length(trackRange)
    
  %%
  currFrameNum = trackRange(i);

  % i is the frame index into lincellnum and trackRange. 
  % trackRange(i) gives the actual frame number that corresponds 
  % to a given lincellnum{i}.
  % Also, recall that a schnitz' frames array is 1-based per 
  % makeschnitz hack to conform to schnitzedit expectations
  
  % empty previous data
  eval([ 'clear Lc ' reg1 ' ' reg2 ' ' reg3 ';']);
  timestamp = 0;

  
  % empty previous data
  %clear Lc yreg rreg;
  %timestamp = 0;

  % load segmented image for this frameNum
  name = [p.segmentationDir,p.movieName,'seg',str3(currFrameNum)];
  disp([' Processing frame ', str3(currFrameNum), ' as nr ', str3(i), ' in range']);
  load([name]); % including variables: LNsub, Lc, phsub, timestamp, rect, and fluor-info, e.g. yback, yreg, yshift, expty, gainy, rback, rreg, rshift, exptr, gainr    

  % get shape properties for cells in this frame
  rp = regionprops(Lc,'MajorAxisLength','MinorAxisLength','Centroid','Orientation','Solidity');

  %------------------------------------------------------------------------
  % LOOP OVER EACH SCHNITZ THAT EXISTS DURING THIS FRAME, AND UPDATE IT
  %------------------------------------------------------------------------
  %%
  schnitzesForFrame = lincellnum{i}; % [schnitznum for cellno1, schnitznum for cellno2, etc]
  nonZeroSchnitzes = schnitzesForFrame(schnitzesForFrame~=0); % with correction you sometimes end up with unexisting schnitzes (0)
  for s = nonZeroSchnitzes
    % figure out index within this schnitz' age-based arrays
    age = find((schnitzcells(s).frame_nrs) == currFrameNum); % MW 2014/06/11 removal N+1 bug
    if isempty(age)
      error(['lincellnum says schnitz num ' s 'exists in frame ' currFrameNum ', but that frame can''t be found in the schnitz' frames array']);
    end

    % get cell number in segmented image for current schnitz & frame
    cellnum = schnitzcells(s).cellno(age);

    % [y,x] are pixels in Lc image where this cell is located
    [y,x] = find(Lc == cellnum); % note: returns (row, column), which will be used as (y,x)
    
    %--------------------------------------------------------------------
    % Add segmentation data: timestamp, areaPixels, area, rp_length, 
    % rp_width, rp_cenX_crop, rp_cenY_crop, rp_cenX_full, rp_cenY_full, 
    % rp_angle, rp_solidity, rp_volume
    %--------------------------------------------------------------------
    schnitzcells(s).timestamp(age)    = timestamp;
    schnitzcells(s).areaPixels(age)   = length(x);
    schnitzcells(s).area(age)         = p.micronsPerPixel * p.micronsPerPixel * length(x);

    % From region properties
    schnitzcells(s).rp_length(age)    = p.micronsPerPixel * rp(cellnum).MajorAxisLength;
    schnitzcells(s).rp_width(age)     = p.micronsPerPixel * rp(cellnum).MinorAxisLength;
    schnitzcells(s).rp_cenX_crop(age) = p.micronsPerPixel * (rp(cellnum).Centroid(1) + rect(2) - 1); % -1 omdat als helemaal links, dan rect(2) = 1
    schnitzcells(s).rp_cenY_crop(age) = p.micronsPerPixel * (rp(cellnum).Centroid(2) + rect(1) - 1); % -1 omdat als helemaal top, dan rect(1) = 1   
    schnitzcells(s).rp_cenX_full(age) = p.micronsPerPixel * (rp(cellnum).Centroid(1) + rect(2) + cropLeftTop(2) - 2); % -2 since twice crop area
    schnitzcells(s).rp_cenY_full(age) = p.micronsPerPixel * (rp(cellnum).Centroid(1) + rect(1) + cropLeftTop(1) - 2); % -2 since twice crop area
    schnitzcells(s).rp_angle(age)     = rp(cellnum).Orientation;
    schnitzcells(s).rp_solidity(age)  = rp(cellnum).Solidity;
                                    r = 0.5 * schnitzcells(s).rp_width(age); % cigar shape estimation of volume of cell in micron^3
    schnitzcells(s).rp_volume(age)    = (4/3)*pi*r*r*r + pi*r*r*(schnitzcells(s).rp_length(age)-2*r);

    
    %--------------------------------------------------------------------
    % Add fitted line: is done, after pixels are rotated so that cell lies straight
    % Using a 2nd and a 3rd degree polynomial
    % Adding: phi, fitCoef2, fitCoef3
    %--------------------------------------------------------------------    
    % rotate pixels, so that cell lies straight
    phi = schnitzcells(s).rp_angle(age)*pi/180; % convert orientation to radians, and use opposite to rotate back
    schnitzcells(s).phi(age) = phi;
    x_rot = x*cos(phi) - y*sin(phi); % mathematical rotation
    y_rot = x*sin(phi) + y*cos(phi); % mathematical rotation

    % Fit 2nd degree polynomial to rotated pixels
    fitCoef2 = DJK_polyfit(x_rot,y_rot,2);
    schnitzcells(s).fitCoef2(age,:) = fitCoef2;

    % Fit 3rd degree polynomial to rotated pixels
    fitCoef3 = DJK_polyfit(x_rot,y_rot,3);
    schnitzcells(s).fitCoef3(age,:) = fitCoef3;
    
     % Fit 4rd degree polynomial to rotated pixels (ADDED)
    fitCoef4 = DJK_polyfit(x_rot,y_rot,4);
    schnitzcells(s).fitCoef4(age,:) = fitCoef4;
    
     % Fit 5rd degree polynomial to rotated pixels (ADDED)
    fitCoef5 = DJK_polyfit(x_rot,y_rot,5);
    schnitzcells(s).fitCoef5(age,:) = fitCoef5;
    
      % Fit 6rd degree polynomial to rotated pixels (ADDED)
    fitCoef6 = DJK_polyfit(x_rot,y_rot,6);
    schnitzcells(s).fitCoef6(age,:) = fitCoef6;
    
      % Fit 7rd degree polynomial to rotated pixels (ADDED)
    fitCoef7 = DJK_polyfit(x_rot,y_rot,7);
    schnitzcells(s).fitCoef7(age,:) = fitCoef7;
    
    func_3rd = @(x)x.^7 .* fitCoef7(1)+x.^6 .* fitCoef7(2)+x.^5 .* fitCoef7(3) + x.^4 .* fitCoef7(4) +x.^3 .* fitCoef7(5) + x.^2 .* fitCoef7(6) + x .* fitCoef7(7) + fitCoef7(8);

    %--------------------------------------------------------------------
    % Determine length from 3rd degree polynomial, OLD WAY
    % Adding: fitCoef3_x_rot_left, fitCoef3_x_rot_right, length_fitCoef3
    %--------------------------------------------------------------------    
    y_line_rot = func_3rd(x_rot);

    % Check which y_rot is closer than 0.5 away from line
    idx_pixel_intersects = find( y_rot>= y_line_rot-0.5 & y_rot<= y_line_rot+0.5 );
    x_rot_intersect = x_rot(idx_pixel_intersects);
    y_rot_intersect = y_rot(idx_pixel_intersects);

    % Determine length
    x_rot_intersect_left  = min(x_rot_intersect);
    x_rot_intersect_right = max(x_rot_intersect);
    func_length = @(x) sqrt( abs( 3 .* x.^2 .* fitCoef3(1) + 2 .* x .* fitCoef3(2) + fitCoef3(3) + 1 ) );
    schnitzcells(s).fitCoef3_x_rot_left(age)   = x_rot_intersect_left;
    schnitzcells(s).fitCoef3_x_rot_right(age)  = x_rot_intersect_right;
    schnitzcells(s).length_fitCoef3(age) = quad(func_length, x_rot_intersect_left,x_rot_intersect_right) * p.micronsPerPixel;

    % Add fluor1_all_frames if color existent
    if strcmp(p.fluor1,'none')==0
         eval(['schnitzcells(s).' fluor1_frames_all '(age)   = NaN;']);
         eval(['existfluor=exist(''' reg1 ''');']) % check if fluorescence picture
         if existfluor==1
             eval(['schnitzcells(s).' fluor1_frames_all '(age) = schnitzcells(s).frame_nrs(age);']); % MW 2014/06/11 removal N+1 bug
             fluor1counter=fluor1counter+1;
         end
    end 
    % Add fluor2_all_frames if color existent
    if strcmp(p.fluor2,'none')==0
         eval(['schnitzcells(s).' fluor2_frames_all '(age)   = NaN;']);
         eval(['existfluor=exist(''' reg2 ''');']) % check if fluorescence picture
         if existfluor==1
             eval(['schnitzcells(s).' fluor2_frames_all '(age) = schnitzcells(s).frame_nrs(age);']); % MW 2014/06/11 removal N+1 bug
             fluor2counter=fluor2counter+1;
         end
    end 
    % Add fluor3_all_frames if color existent
    if strcmp(p.fluor3,'none')==0
         eval(['schnitzcells(s).' fluor3_frames_all '(age)   = NaN;']);
         eval(['existfluor=exist(''' reg3 ''');']) % check if fluorescence picture
         if existfluor==1
             eval(['schnitzcells(s).' fluor3_frames_all '(age) = schnitzcells(s).frame_nrs(age);']); % MW 2014/06/11 removal N+1 bug
             fluor3counter=fluor3counter+1;
         end
    end 
    
  end % loop over schnitzesForFrame
end % loop over frames
%--------------------------------------------------------------------------

  
%--------------------------------------------------------------------------
% LOOP OVER SCHNITZCELLS TO ADD FINAL STUFF
%--------------------------------------------------------------------------
% determine timestamp for time=0
all_timestamps = [schnitzcells.timestamp];
s = sort([all_timestamps(find(all_timestamps>0))]);
firstmoment = s(1);

% loop over schnitzcells
for i = 1:length(schnitzcells)
    
   %--------------------------------------------------------------------
  % Add time data: time, av_time,  
  % completeCycle, birthTime, divTime, interDivTime, gen, phase
  %--------------------------------------------------------------------
  % add time
  tdiff = schnitzcells(i).timestamp - firstmoment;
  schnitzcells(i).time    = DJK_getMinutesFromTimestamp(tdiff);
  schnitzcells(i).av_time = schnitzcells(i).time(1) + 0.5*(schnitzcells(i).time(end) - schnitzcells(i).time(1));
   
  %--------------------------------------------------------------------
  % Convert fluor1_frames_all to fluor_frames  (e.g. Y_frames_all to
  % Y_frames) (if color existent)
  % Add time for fluor frames: fluor1_time. fluor2_time, fluor3_time (e.g. Y_time),
  %--------------------------------------------------------------------
  if (strcmp(p.fluor1,'none')==0)
    eval(['idx = find(~isnan(schnitzcells(i).' fluor1_frames_all '));']);
    eval(['schnitzcells(i).' fluor1_frames '= schnitzcells(i).' fluor1_frames_all '( idx );']);
    eval(['schnitzcells(i).' fluor1_time ' = schnitzcells(i).time( idx );']);
  end
  %--------------------------------------------------------------------
  % Convert fluor2_frames_all to fluor2_frames
  %--------------------------------------------------------------------
  if (strcmp(p.fluor2,'none')==0)
    eval(['idx = find(~isnan(schnitzcells(i).' fluor2_frames_all '));']);
    eval(['schnitzcells(i).' fluor2_frames '= schnitzcells(i).' fluor2_frames_all '( idx );']);
    eval(['schnitzcells(i).' fluor2_time ' = schnitzcells(i).time( idx );']);
  end
  %--------------------------------------------------------------------
  % Convert fluor3_frames_all to fluor3_frames
  %--------------------------------------------------------------------
  if (strcmp(p.fluor3,'none')==0)
    eval(['idx = find(~isnan(schnitzcells(i).' fluor3_frames_all '));']);
    eval(['schnitzcells(i).' fluor3_frames '= schnitzcells(i).' fluor3_frames_all '( idx );']);
    eval(['schnitzcells(i).' fluor3_time ' = schnitzcells(i).time( idx );']);
  end
  %--------------------------------------------------------------------
  
  
  % A cell has a complete cell cycle when it's birth and division is seen
  schnitzcells(i).completeCycle = 0;
  if (schnitzcells(i).P ~= 0 & schnitzcells(i).E ~= 0 & schnitzcells(i).D ~= 0 )
    schnitzcells(i).completeCycle = 1;
  end

  % Birth, division and interdivision time only in case of available data
  P = schnitzcells(i).P; E = schnitzcells(i).E; D = schnitzcells(i).D;
  schnitzcells(i).birthTime = NaN; schnitzcells(i).divTime = NaN; schnitzcells(i).interDivTime = NaN;
  if P
    schnitzcells(i).birthTime = DJK_getMinutesFromTimestamp(0.5*(schnitzcells(i).timestamp(1) + schnitzcells(P).timestamp(end)) - firstmoment);
%     schnitzcells(i).birthTime = 0.5*(schnitzcells(i).time(1) + schnitzcells(P).time(end));
  end
  if E & D
    schnitzcells(i).divTime   = DJK_getMinutesFromTimestamp(0.5*(schnitzcells(i).timestamp(end) + schnitzcells(E).timestamp(1)) - firstmoment);
%     schnitzcells(i).divTime   = 0.5*(schnitzcells(i).time(end) + schnitzcells(E).time(1));
  end
  if (schnitzcells(i).completeCycle)
%     schnitzcells(i).interDivTime = DJK_getMinutesFromTimestamp(schnitzcells(i).divTime - schnitzcells(i).birthTime);    
    schnitzcells(i).interDivTime = schnitzcells(i).divTime - schnitzcells(i).birthTime;    
  end

  % Add generations
  gen = 0; me = i;
  while (schnitzcells(me).P ~= 0)
    gen = gen+1;
    me = schnitzcells(me).P;
  end
  schnitzcells(i).gen = gen;
  
  % Add phase
  schnitzcells(i).phase = NaN * schnitzcells(i).frame_nrs;
  if (schnitzcells(i).completeCycle)
    schnitzcells(i).phase = (schnitzcells(i).time - schnitzcells(i).birthTime) / schnitzcells(i).interDivTime;
  end

  % Add sister
  schnitzcells(i).S = 0; 
  if P
    if i ~= schnitzcells(P).D, schnitzcells(i).S = schnitzcells(P).D; end
    if i ~= schnitzcells(P).E, schnitzcells(i).S = schnitzcells(P).E; end
  end    
end
%--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % REMOVE FIELDS I WILL PROBABLY NOT USE
% %--------------------------------------------------------------------------
% rmfields = {''};
% schnitzcells = rmfield(schnitzcells, rmfields);
% %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Display some information
%--------------------------------------------------------------------------
% Count nr of cells with complete cell cycle
completeCycleNr = 0;

% loop over schnitzcells
for i = 1:length(schnitzcells)
  if schnitzcells(i).completeCycle
    completeCycleNr = completeCycleNr+1;
  end
end

disp(['-----------------------------------------------------------']);
% nr of cells with complete cell cycle
disp([' * Total nr cells with complete cell cycle in schnitzcells: ', num2str(completeCycleNr)]);
disp(['-----------------------------------------------------------']);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Display warning if no fluorescence pictures were found
%--------------------------------------------------------------------------
if (strcmp(p.fluor1,'none')==0 & fluor1counter==0)
    disp(['Warning! No fluorescence pictures for fluor1 (' p.fluor1 ') found. ' ...
          'Maybe color in ''p'' wrong associated. Change with p.fluor1=...']);
end
if (strcmp(p.fluor2,'none')==0 & fluor2counter==0)
    disp(['Warning! No fluorescence pictures for fluor2 (' p.fluor2 ') found. ' ...
          'Maybe color in ''p'' wrong associated. Change with p.fluor2=...']);
end
if (strcmp(p.fluor3,'none')==0 & fluor3counter==0)
    disp(['Warning! No fluorescence pictures for fluor3 (' p.fluor3 ') found. ' ...
          'Maybe color in ''p'' wrong associated. Change with p.fluor3=...']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Save extacted data
%--------------------------------------------------------------------------
save(p.schnitzName, 'schnitzcells');
disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------
