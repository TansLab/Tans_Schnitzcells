function [p,schnitzcells] = DJK_addToSchnitzes_fluor_anycolor(p,varargin) 
% DJK_addToSchnitzes_fluor_anycolor loads schnitzcells, calculates the fluor for
% schnitzes and saves them to schnitzcells again.
%
% adjust name of 'autoFluor' to be color specific??? Not until now.
%
% DJK_compileSchnitzImproved & DJK_addToSchnitzes_length must have been run before 
%
% OUTPUT
% 'p'   
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'p'
%
% OPTIONAL ARGUMENTS:
% 'fluorcolor'    Fluor color that will be added. This can be
%                 'fluor1', 'fluor2', 'fluor3'. The colors associated with
%                 these expressions are stored in 'p'. If no color is
%                 given, 'fluor1'will be chosen.
% 'schnitzName'       
% 'autoFluor2'         auto fluorescence of cells, subtracted from fluor
%                     default: = 0
% 'micronsPerPixel'   default = 0.04065 (1.5x magnification of Killer Mike)
% 'blockSettings'     used for blockFluor determination: indicates distance
%                     from poles, width of block, & max distance offset, all in um
%                     default: [0.325 0.4 0.06]
% 'onScreen' = 0      automatically save and close images
%              1      will ask to save and close images
%              2      will not show or save images (default)
% 'colorNormalize'    default: [0 15]
% 'DJK_saveDir2'       Directory where images should be saved. Defaults to "p.analysisDir \ fluor2 \ R6fromBlock \"
% 'minimalMode'       =1: only calculates values, that are actually used
%                     later (Y5xxx,Y6_mean,(+ElowitzStyle) etc).
%                     =0: calculates everything (default: =0)  (NW 2012/04)
%
% NOTES ON PROCEDURE AND OUTPUT
% ===
% In principle, pixels are selected from the fluor image using the 
% segmentation file, which contains the regions of the detected cells
% encoded as indices in a matrix. 
% There are however different normalizations (e.g. illuimation time)
% and corrections performed on the fluor images,
% and one can determine either the sum or the average of the fluor
% intensity within a cell. Furthermore, an additional method is introduced
% here to take only the fluor values from a subset of the detected cell
% area, namely only its central area, to prevent artifacts introduced by
% (incorrect) border detection.
% All these options lead to fluor values being calculated differently, and
% these different values are all stored.
% The parameters <fluor>1_suffix to <fluor>5_suffix (where <fluor> is a 
% capital letter encoding for the fluor, e.g. G, C, Y, ..) hold the 
% different consecutive corrections that were performed on the fluor image. 
% fluor6 holds the fluor value determined from the central area. The
% suffixes of these parameters tell you by which method the different pixel
% intensities are summarized (e.g. sum, mean, ..).

DEBUG = 1;

%--------------------------------------------------------------------------
%% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_addToSchnitzes_fluor_anycolor';

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
    % fluorcolor is special case, because the CONTENT of p.fluor1/2/3 must be written to
    % p.fluorcolor and not fluor1/2/3 literally. (NW 11/12/08)
    if strcmp(fieldName,'fluorcolor')==1 
        p.(fieldName)=p.(varargin{i+1});
        disp('--------------------------------------------------------------');
        disp(['Using ' varargin{i+1} ' (' p.(fieldName) ') as fluorescence color.']);
    else
        p.(fieldName) = varargin{i+1};
    end;
  end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% Set default parameter values if they don't exist yet
if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end
if ~existfield(p,'autoFluor2')
  p.autoFluor2 = 0;
end
if ~existfield(p,'micronsPerPixel')
  p.micronsPerPixel = 0.04065;
end
if ~existfield(p,'blockSettings')
  p.blockSettings = [0.325 0.2 0.06]; % [0.325 0.2 0.06];
end
% convert p.blockSettings to valuable pixel values
pole_dist   = p.blockSettings(1) / p.micronsPerPixel;
block_width = p.blockSettings(2) / p.micronsPerPixel;
max_offset  = p.blockSettings(3) / p.micronsPerPixel;

if ~existfield(p,'onScreen')
  p.onScreen = 2;
end
if ~existfield(p,'colorNormalize')
  p.colorNormalize = [0 15];
end
if ~existfield(p,'minimalMode') % (NW 2012/04)
  p.minimalMode = 0;
end

% Choose fluorescence color if not specified
if ~isfield(p,'fluorcolor')
    p.fluorcolor=p.fluor1;
    disp('--------------------------------------------------------------');
    disp(['No fluorescence color specified. Will take fluor1 (' p.fluor1 ').']);
end

% If explicit DJK_saveDir2 is not given, define it
if ~existfield(p,'DJK_saveDir2')
  p.DJK_saveDir2 = [p.analysisDir 'fluor_' p.fluorcolor filesep upper(p.fluorcolor) 'fromFluor' filesep];
end
% make sure every directory field has a trailing filesep
pathlength = length(p.DJK_saveDir2);
if (p.DJK_saveDir2(pathlength) ~= filesep)
  p.DJK_saveDir2 = [p.DJK_saveDir2 filesep];
end
% Make sure DJK_saveDir2 exists, else make it
if exist(p.DJK_saveDir2)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir2]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir2 ' : ' msg]);
    return;
  end
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Load Schnitzcells
%--------------------------------------------------------------------------
% Loading an existing schnitz file which contains the image-derived fields
if exist(p.schnitzName)~=2
  error(['Could not read ' p.schnitzName ' , which is required for quick mode']);
end
load(p.schnitzName);
disp(['Load from ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% CREATION OF lincellnum
%--------------------------------------------------------------------------
% Will use a lincellnum structure, where for only used frames (not necessarly
% starting with image 001 (frame 2), the link between schnitznum and cellno
% is saved: lincellnum {a} (b) returns schnitznum of cellno b in frame a

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
  for age = 1:length(s.frame_nrs)
    framenum = s.frame_nrs(age); % MW 2014/06/11 removal N+1 bug

    % look up that frame number's index within trackRange
    lincellnumIndex = find(trackRange==framenum);
    cellnum  = s.cellno(age);
    lincellnum{lincellnumIndex}(cellnum) = schnitznum;
  end
end
%--------------------------------------------------------------------------


%%--------------------------------------------------------------------------
% check if fluorcolor is non-existent ('none')
%--------------------------------------------------------------------------
if strcmp(p.fluorcolor,'none')==1
    disp('Fluorescence color is non-existent (set to ''none'').');
    %optimalShift=[];
    return
end

%%

% -------------------------------------------------
% create variable names
% -------------------------------------------------
reg=genvarname([p.fluorcolor 'reg']);         % e.g. yreg
binning=genvarname([p.fluorcolor 'binning']); % e.g. ybinning
back=genvarname([p.fluorcolor 'back']);       % e.g. yback
gain=genvarname(['gain' p.fluorcolor]);       % e.g. gainy
expt=genvarname(['expt' p.fluorcolor]);       % e.g. expty
reg2=genvarname([p.fluorcolor 'reg2']);           % e.g. yreg2   %HERE WAS ACCIDENTLY '5' IN ALL VARIABLES!!!! (NW 2012-04-26)
reg3=genvarname([p.fluorcolor 'reg3']);           % e.g. yreg3    %HERE WAS ACCIDENTLY '5' IN ALL VARIABLES!!!!
reg4=genvarname([p.fluorcolor 'reg4']);           % e.g. yreg4    %HERE WAS ACCIDENTLY '5' IN ALL VARIABLES!!!!
reg5=genvarname([p.fluorcolor 'reg5']);           % e.g. yreg5 
regNorm=genvarname([p.fluorcolor 'regNorm']);     % e.g. yregNorm 
reg2Norm=genvarname([p.fluorcolor 'reg2Norm']);   % e.g. yreg2Norm 
reg3Norm=genvarname([p.fluorcolor 'reg3Norm']);   % e.g. yreg3Norm 
reg4Norm=genvarname([p.fluorcolor 'reg4Norm']);   % e.g. yreg4Norm 
reg5Norm=genvarname([p.fluorcolor 'reg5Norm']);   % e.g. yreg5Norm 
back2=genvarname([p.fluorcolor 'back2']);
back3=genvarname([p.fluorcolor 'back3']);
back2Alt=genvarname([p.fluorcolor 'back2Alt']);
fluor_frames_all=genvarname([upper(p.fluorcolor) '_frames_all']);
fluor_sum_all=genvarname([upper(p.fluorcolor) '_sum_all']);
fluor_mean_all=genvarname([upper(p.fluorcolor) '_mean_all']);
fluor_back_all=genvarname([upper(p.fluorcolor) '_back_all']);
fluor2_sum_all=genvarname([upper(p.fluorcolor) '2_sum_all']);
fluor2_mean_all=genvarname([upper(p.fluorcolor) '2_mean_all']);
fluor2_back_all=genvarname([upper(p.fluorcolor) '2_back_all']);
fluor2_backAlt_all=genvarname([upper(p.fluorcolor) '2_backAlt_all']);
fluor3_sum_all=genvarname([upper(p.fluorcolor) '3_sum_all']);
fluor3_mean_all=genvarname([upper(p.fluorcolor) '3_mean_all']);
fluor3_back_all=genvarname([upper(p.fluorcolor) '3_back_all']);
fluor4_sum_all=genvarname([upper(p.fluorcolor) '4_sum_all']);
fluor4_mean_all=genvarname([upper(p.fluorcolor) '4_mean_all']);
fluor5_sum_all=genvarname([upper(p.fluorcolor) '5_sum_all']);
fluor5_mean_all=genvarname([upper(p.fluorcolor) '5_mean_all']);
fluor5_stdev_all=genvarname([upper(p.fluorcolor) '5_stdev_all']);
fluor5_median_all=genvarname([upper(p.fluorcolor) '5_median_all']);
fluor_sum=genvarname([upper(p.fluorcolor) '_sum']);
fluor_mean=genvarname([upper(p.fluorcolor) '_mean']);
fluor_back=genvarname([upper(p.fluorcolor) '_back']);
fluor2_sum=genvarname([upper(p.fluorcolor) '2_sum']);
fluor2_mean=genvarname([upper(p.fluorcolor) '2_mean']);
fluor2_back=genvarname([upper(p.fluorcolor) '2_back']);
fluor2_backAlt=genvarname([upper(p.fluorcolor) '2_backAlt']);
fluor3_sum=genvarname([upper(p.fluorcolor) '3_sum']);
fluor3_mean=genvarname([upper(p.fluorcolor) '3_mean']);
fluor3_back=genvarname([upper(p.fluorcolor) '3_back']);
fluor4_sum=genvarname([upper(p.fluorcolor) '4_sum']);
fluor4_mean=genvarname([upper(p.fluorcolor) '4_mean']);
fluor5_sum=genvarname([upper(p.fluorcolor) '5_sum']);
fluor5_mean=genvarname([upper(p.fluorcolor) '5_mean']);
fluor5_stdev=genvarname([upper(p.fluorcolor) '5_stdev']);
fluor5_median=genvarname([upper(p.fluorcolor) '5_median']);
fluor6_mean_all=genvarname([upper(p.fluorcolor) '6_mean_all']);
fluor6_offset_all=genvarname([upper(p.fluorcolor) '6_offset_all']);
fluor6_mean=genvarname([upper(p.fluorcolor) '6_mean']);
fluor_time=genvarname([upper(p.fluorcolor) '_time']);
temp_fluor_mean=genvarname(['temp_' upper(p.fluorcolor) '_mean']);
temp_fitted_fluor_mean=genvarname(['temp_fitted_' upper(p.fluorcolor) '_mean']);
% -------------------------------------------------

% --------------------------------------------------
% create control variables for check if right fluorescence colors were
% added to 'p'
% --------------------------------------------------
fluorcounter=0; 
% --------------------------------------------------

%--------------------------------------------------------------------------
%% LOOP OVER FRAMES IN TRACKRANGE
%--------------------------------------------------------------------------
for num = 1:length(trackRange)
    
  %%
  fr = trackRange(num);

  %%
  % empty previous data
  eval(['clear Lc ' reg ' ' reg4 ' ' reg5 ' ' reg5Norm ';']);

  % load segmented image for this frameNum
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load([name]); % including variables: LNsub, Lc, phsub, timestamp, rect, and e.g. yback, yreg, yshift, expty, gainy, rback, rreg, rshift, exptr, gainr    

  % in case onScreen, Lc_block will be used to save where blocks are
  Lc_block = zeros(size(Lc));

  %------------------------------------------------------------------------
  % CHECK WHETHER FLUOR EXIST FOR THIS FRAME
  %------------------------------------------------------------------------
  eval(['existfluor=exist(''' reg ''');']) % check if fluorescence picture exists (=1)
  if existfluor==0
     disp([' * Frame ' str3(fr) ' : no fluor (' p.fluorcolor ') image']);
    % skip frame if no fluor
  else
    fluorcounter=fluorcounter+1;

    disp([' * Frame ' str3(fr) ' : processing as nr ' str3(num) ' in range']);

    % load fluor for this frameNum
    name= [p.tracksDir, p.movieName, 'Fluor_', p.fluorcolor, '_' , str3(fr)];
    disp(['   -> loading Fluor of ', str3(fr)]);
    load([name]); % including variables: exptr, gainr, phaseCropSize, phaseFullSize, rect, rectCrop, and e.g. rback, rback2, rback2Alt, rbackAlt, rbinning, rreg, rreg2

    %----------------------------------------------------------------------
    % NORMALIZE FLUOR IMAGES
    %----------------------------------------------------------------------
    % reg    : gedeelte van fluor image, vergroot naar gelang binning
    % reg2   : flatfield & shading corrected reg
    % reg3   : shift corrected reg2
    % reg4   : convolved version of reg2
    % reg5   : shift corrected reg4
    % back   : background of reg, traditional
    % back2  : background of reg2, traditional
    % back3  : background of reg3, traditional

    % normalize fluor images: correct for background, exposure time, binning, micronsPerPixel & autoFluor2
    eval([ regNorm '  = normalizeElowitz( ' reg ',' back ',' expt ',' binning ',p.autoFluor2);']);

    if p.minimalMode~=1
        eval([ reg2Norm  ' = normalize(' reg2 ',' back2 ',' expt ',' binning ', p.micronsPerPixel,p.autoFluor2);']);
        eval([ reg3Norm ' = normalize(' reg3 ',' back3 ',' expt ',' binning ', p.micronsPerPixel,p.autoFluor2);']);
    end

    % if deconvolved image exist do the same for them
    eval(['existfluor=exist(''' reg5 ''');']) % check if deconvolved picture exists (=1) %take reg5 instead of reg4 to account for minimalMode (NW 2012-04, doesn't calculate everything)
    if (existfluor==1)
        if p.minimalMode~=1
            eval([ reg4Norm '= normalize(' reg4 ',' back2 ',' expt ',' binning ',p.micronsPerPixel,p.autoFluor2);']);
        end
        eval([ reg5Norm '= normalize(' reg5 ',' back3 ',' expt ',' binning ',p.micronsPerPixel,p.autoFluor2);']);
    end
    %----------------------------------------------------------------------
  end
  %------------------------------------------------------------------------

  %%
  %----------------------------------------------------------------------
  % LOOP OVER EACH SCHNITZ THAT EXISTS DURING THIS FRAME, AND UPDATE IT
  %----------------------------------------------------------------------
  schnitzesForFrame = lincellnum{num}; % [schnitznum for cellno1, schnitznum for cellno2, etc]
  nonZeroSchnitzes = schnitzesForFrame(schnitzesForFrame~=0); % with correction you sometimes end up with unexisting schnitzes (0)

  for i = nonZeroSchnitzes
    % figure out index within this schnitz' age-based arrays
    age     = find((schnitzcells(i).frame_nrs) == fr); % MW 2014/06/11 removal N+1 bug
    cellno  = schnitzcells(i).cellno(age);

    if isempty(age)
      error(['lincellnum says schnitz num ' i 'exists in frame ' fr ', but that frame can''t be found in the schnitz' frames array']);
    end

    schnitzcells(i).autoFluor2_all(age)     = p.autoFluor2;

    %--------------------------------------------------------------------
    % First fill fluor data with NaN
    %--------------------------------------------------------------------
    eval(['schnitzcells(i).' fluor_sum_all '(age)        = NaN;'])
    eval(['schnitzcells(i).' fluor_mean_all '(age)       = NaN;'])
    eval(['schnitzcells(i).' fluor_back_all '(age)       = NaN;'])

    if p.minimalMode~=1
        eval(['schnitzcells(i).' fluor2_sum_all '(age)       = NaN;'])
        eval(['schnitzcells(i).' fluor2_mean_all '(age)      = NaN;'])
        eval(['schnitzcells(i).' fluor2_back_all '(age)      = NaN;'])
        eval(['schnitzcells(i).' fluor2_backAlt_all '(age)   = NaN;'])

        eval(['schnitzcells(i).' fluor3_sum_all '(age)       = NaN;'])
        eval(['schnitzcells(i).' fluor3_mean_all '(age)      = NaN;'])
        eval(['schnitzcells(i).' fluor3_back_all '(age)      = NaN;'])

        eval(['schnitzcells(i).' fluor4_sum_all '(age)       = NaN;'])
        eval(['schnitzcells(i).' fluor4_mean_all '(age)      = NaN;'])
    end

    eval(['schnitzcells(i).' fluor5_sum_all '(age)       = NaN;'])
    eval(['schnitzcells(i).' fluor5_mean_all '(age)      = NaN;'])
    eval(['schnitzcells(i).' fluor5_stdev_all '(age)     = NaN;'])
    eval(['schnitzcells(i).' fluor5_median_all '(age)    = NaN;'])
    %--------------------------------------------------------------------


    %--------------------------------------------------------------------
    % Add fluor data:
    % R_sum_all, SR, R_mean_all, medR, backR, Rall_sum, Rall_stdev, Rall_mean, 
    % Rall_median, Rall_back, Rall_backAlt
    %--------------------------------------------------------------------    
    eval(['existfluor=exist(''' reg ''');']) % check if fluorescence picture exists (=1) 
    if existfluor==1

      % get cell number in segmented image for current schnitz & frame
      cellnum = schnitzcells(i).cellno(age);

      % [loc_x, loc_y] are pixels in Lc image where this cell is located
      loc = find(Lc == cellnum); 

      % calc Traditional Fluor data
      eval(['schnitzcells(i).' fluor_sum_all '(age)        = sum(' regNorm '(loc));'])
      eval(['schnitzcells(i).' fluor_mean_all '(age)       = schnitzcells(i).' fluor_sum_all '(age)/schnitzcells(i).areaPixels(age);'])
      eval(['schnitzcells(i).' fluor_back_all '(age)       = ' back ';'])

      if p.minimalMode~=1
          % calc Flatfield/Shading corrected Fluor data
          eval(['schnitzcells(i).' fluor2_sum_all '(age)       = sum(' reg2Norm '(loc));'])
          eval(['schnitzcells(i).' fluor2_mean_all '(age)      = schnitzcells(i).' fluor2_sum_all '(age)/schnitzcells(i).areaPixels(age);'])
          eval(['schnitzcells(i).' fluor2_back_all '(age)      = ' back2 ';'])
          eval(['schnitzcells(i).' fluor2_backAlt_all '(age)   = ' back2Alt ';'])

          % calc Flatfield/Shading/Shift corrected Fluor data
          eval(['schnitzcells(i).' fluor3_sum_all '(age)       = sum(' reg3Norm '(loc));'])
          eval(['schnitzcells(i).' fluor3_mean_all '(age)      = schnitzcells(i).' fluor3_sum_all '(age)/schnitzcells(i).areaPixels(age);'])
          eval(['schnitzcells(i).' fluor3_back_all '(age)      = ' back3 ';'])
      end

      eval(['existfluor=exist(''' reg5 ''');']) % check if fluorescence picture exists (=1) %take reg5 instead of reg4 to account for minimalMode (NW 2012-04, doesn't calculate everything)
      if existfluor==1
        % calc Flatfield/Shading/Deconvolved corrected Fluor data
        if p.minimalMode~=1
            eval(['schnitzcells(i).' fluor4_sum_all '(age)     = sum(' reg4Norm '(loc));'])
            eval(['schnitzcells(i).' fluor4_mean_all '(age)    = schnitzcells(i).' fluor4_sum_all '(age)/schnitzcells(i).areaPixels(age);'])
        end

        eval(['schnitzcells(i).' fluor5_sum_all '(age)     = sum(' reg5Norm '(loc));'])
        eval(['schnitzcells(i).' fluor5_mean_all '(age)    = schnitzcells(i).' fluor5_sum_all '(age)/schnitzcells(i).areaPixels(age);'])
        eval(['schnitzcells(i).' fluor5_stdev_all '(age)   = std(' reg5Norm '(loc));'])
        eval(['schnitzcells(i).' fluor5_median_all '(age)  = median(' reg5Norm '(loc));'])
      end
    end


    %----------------------------------------------------------------------
    % Add block fluor of rreg5Norm
    %----------------------------------------------------------------------
    eval(['schnitzcells(i).' fluor6_mean_all '(age)      = NaN;'])
    eval(['schnitzcells(i).' fluor6_offset_all '(age)    = NaN;'])

    eval(['existfluor=exist(''' reg5Norm ''');']) % check if variable exists (=1)
    if existfluor==1
      % get rotated pixels, so that cell lies straight
      [y,x] = find(Lc==cellno); % note: returns (row, column), which will be used as (y,x)
      phi = schnitzcells(i).rp_angle(age)*pi/180; % convert orientation to radians, and use opposite to rotate back
      x_rot = x*cos(phi) - y*sin(phi); % mathematical rotation
      y_rot = x*sin(phi) + y*cos(phi); % mathematical rotation

      % Get 3rd degree polynomial to rotated pixels
      fitCoef3            = schnitzcells(i).fitCoef3(age,:);
      fitNew_x_rot_left   = schnitzcells(i).fitNew_x_rot_left(age);
      fitNew_x_rot_right  = schnitzcells(i).fitNew_x_rot_right(age);

      offset_measurements = struct; % will contain data

      % determine offset values
      temp = 0.5*floor(max_offset/0.5);
      offset_values = [-max_offset -temp:0.5:temp max_offset];

      % loop over offset values
      for offset_nr = 1:length(offset_values)
        % offsetted 
        func_3rd = @(x) x.^3 .* fitCoef3(1) + x.^2 .* fitCoef3(2) + x .* fitCoef3(3) + fitCoef3(4) + offset_values(offset_nr);

        % determine points on fitted line
        x_rot_line = [fitNew_x_rot_left:fitNew_x_rot_right];
        y_rot_line = func_3rd(x_rot_line);

        % Remove pixels that are too close to poles
        fitNew_y_rot_left   = func_3rd(fitNew_x_rot_left);
        fitNew_y_rot_right  = func_3rd(fitNew_x_rot_right);
        dist_squared_left   = (x_rot_line-fitNew_x_rot_left).^2 + (y_rot_line-fitNew_y_rot_left).^2;
        dist_squared_right  = (x_rot_line-fitNew_x_rot_right).^2 + (y_rot_line-fitNew_y_rot_right).^2;
        min_distance        = (pole_dist+0.5*block_width)^2;
        idx_far_away        = find( dist_squared_left > min_distance & dist_squared_right > min_distance);
        x_rot_line          = x_rot_line(idx_far_away);
        y_rot_line          = y_rot_line(idx_far_away);

        % Determine pixels that are close enough to line
        min_distance        = (0.5*block_width)^2;
        idx_notBlock        = [1:length(x_rot)];
        idx_block           = [];
        for j = 1:length(x_rot_line)
          dist_squared      = (x_rot(idx_notBlock)-x_rot_line(j)).^2 + (y_rot(idx_notBlock)-y_rot_line(j)).^2;
          idx_close         = find( dist_squared <= min_distance);
          idx_block         = [idx_block idx_notBlock(idx_close')];
          idx_notBlock( idx_close' ) = [];
        end

        % Rotate block back
        x_block             = x_rot(idx_block)*cos(-phi) - y_rot(idx_block)*sin(-phi); % mathematical rotation
        y_block             = x_rot(idx_block)*sin(-phi) + y_rot(idx_block)*cos(-phi); % mathematical rotation

        % Save values
        offset_measurements(offset_nr).block = sub2ind(size(Lc), round(y_block), round(x_block)) ;
        eval(['offset_measurements(offset_nr).' fluor6_mean '= mean( ' reg5Norm '( offset_measurements(offset_nr).block ) );'])
%         disp([str3(offset_nr) ' - ' num2str(offset_values(offset_nr)) ' : ' num2str(offset_measurements(offset_nr).R6_mean)]);
      end
      eval(['[trash, best_offset_nr] = max( [offset_measurements(:).' fluor6_mean '] );'])

%       disp(['schnitz ' str3(i) ' frame ' str3(fr) ' best_offset_nr ' str3(best_offset_nr) ' best_offset ' num2str(offset_values(offset_nr)) ' value ' num2str(offset_measurements(offset_nr).R6_mean)]);

      eval(['schnitzcells(i).' fluor6_mean_all '(age)  = offset_measurements(best_offset_nr).' fluor6_mean ';'])
      eval(['schnitzcells(i).' fluor6_offset_all '(age)    = offset_values(best_offset_nr)* p.micronsPerPixel;'])

      %--------------------------------------------------------------------
      % In case onScreen, save location of block
      %--------------------------------------------------------------------
      if p.onScreen ~= 2
        Lc_block( offset_measurements(best_offset_nr).block ) = cellno;
      end
      %--------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
  end % loop over schnitzesForFrame


  %--------------------------------------------------------------------------
  % Make plot
  %--------------------------------------------------------------------------
  eval(['existfluor=exist(''' reg5Norm ''');']) % check if variable exists (=1)

  if p.onScreen ~= 2 & existfluor==1




    %--------------------------------------------------------------------------
    % MAKE COLOR IMAGE OF FLUOR IMAGE
    %--------------------------------------------------------------------------
    % Note: in color_map 1 to 251, corresponds to 0 to 250 in unint16 image
    % Note: in color_map 1 to 251, corresponds to 1 to 251 in double image
    % colormap, will use double image, so image(x,y) = 1, links to color_map(1)
    number_of_colors = 250;
    color_map = hot(number_of_colors); % used to be jet colorGray
    color_text = number_of_colors+1;
    color_map(color_text,:) = [1 1 1]; % for text at colorbar

    % scale fluorimage %blubb
    eval(['fluor = round( DJK_scaleRange(' reg5Norm ', p.colorNormalize, [0 number_of_colors-1]) );'])
    %--------------------------------------------------------------------------



    %--------------------------------------------------------------------------
    % MAKE PERIMETER OF SEG IMAGE
    %--------------------------------------------------------------------------
    Lc_perim = zeros(size(Lc));
    % loop over cells
    for c = 1:max2(Lc)
      Lcell = +(Lc == c);
      Lcell_perim = bwperim(Lcell);
      Lc_perim( find(Lcell_perim>0) ) = c;
    end
    %--------------------------------------------------------------------------

    %----------------------------------------------------------------------
    % MAKE COLOR BAR
    %----------------------------------------------------------------------
    % colorbar on left-top from 1 to 100 (color_map goes from 1 to 100)
    Lc_bar = zeros([110 50]);
    for i=[0:10] % 1 block for 11 values: 100 till 0
      Lc_bar((10*i)+1:(10*(i+1)) , 1:6) = round( DJK_scaleRange(10-i, [0 10], [0 number_of_colors-1]) );
    end
    min_im = DJK_getImageOfNr( p.colorNormalize(1) );
    max_im = DJK_getImageOfNr( p.colorNormalize(2) );
    min_im(find(min_im>0)) = color_text;
    max_im(find(max_im>0)) = color_text;
    Lc_bar = DJK_imagePlace(Lc_bar, min_im, 8+size(min_im,2)/2, 105);
    Lc_bar = DJK_imagePlace(Lc_bar, max_im, 8+size(max_im,2)/2, 5);
    %----------------------------------------------------------------------

    %--------------------------------------------------------------------------
    % ADD PERIMETER & BLOCK & COLOR BAR TO FLUOR
    %--------------------------------------------------------------------------
    fluor_block = fluor;
    fluor_block( find(Lc_perim>0) ) = color_text;
    fluor_block( find(Lc_block>0) ) = 1;  %displays black Y6_block blubb
    fluor_block = DJK_imagePlace(fluor_block, Lc_bar, 25, 55); %adds colorbar  into fluorimage blubb
    %--------------------------------------------------------------------------

    % START. activate, if you want image on black background which has the same
     % size as phase images
%  	dummyimage=ones(phaseCropSize);
%   
%    dummyimage(rect(1):rect(3),rect(2):rect(4))=fluor_block;
%    fluor_block=dummyimage;
     % END


    %--------------------------------------------------------------------
    % Show figure
    %--------------------------------------------------------------------    
    if p.onScreen == 1
      figureName = [p.movieDate ' ' p.movieName ' fluor_' p.fluorcolor ' block fr ' str3(fr)];
      fig11 = figure('color','k', 'Name', figureName);
      imshow(ind2rgb(fluor_block,color_map),'InitialMagnification',100);
    end
    %--------------------------------------------------------------------    

    %----------------------------------------------------------------------
    % ASKING TO SAVE FIGURE
    %----------------------------------------------------------------------
    % Ask to save the figure
    if p.onScreen == 1
      saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
      pause(0.2);
    else
      saveFigInput='Yes and Close';
    end
    if (upper(saveFigInput(1))=='Y')
      figureFileName = [p.movieName '-' p.fluorcolor '6-' str3(fr)];
      imwrite(ind2rgb(fluor_block,color_map),[p.DJK_saveDir2 figureFileName '.png'],'png');
      if (strcmp(saveFigInput,'Yes and Close')) & p.onScreen == 1
        close(fig11);
        pause(0.2);
      end
      disp([' * Saved plot in ' figureFileName '.png']);
    end
    %----------------------------------------------------------------------
  end

end % loop over frames
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZCELLS TO CONVERT _ALL FIELDS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)

  % Find fluor indexes where value exists (i.e. not NaN).
  eval(['FluorIndex = find(~isnan(schnitzcells(i).' fluor_frames_all '));'])

  eval(['schnitzcells(i).autoFluor2       =  schnitzcells(i).autoFluor2_all( FluorIndex );'])

  eval(['schnitzcells(i).' fluor_sum          ' = schnitzcells(i).' fluor_sum_all '( FluorIndex );'])
  eval(['schnitzcells(i).' fluor_mean         ' = schnitzcells(i).' fluor_mean_all '( FluorIndex );'])
  eval(['schnitzcells(i).' fluor_back         ' = schnitzcells(i).' fluor_back_all '( FluorIndex );'])

  if p.minimalMode~=1
      eval(['schnitzcells(i).' fluor2_sum         ' = schnitzcells(i).' fluor2_sum_all '( FluorIndex );'])
      eval(['schnitzcells(i).' fluor2_mean        ' = schnitzcells(i).' fluor2_mean_all '( FluorIndex );'])
      eval(['schnitzcells(i).' fluor2_back        ' = schnitzcells(i).' fluor2_back_all '( FluorIndex );'])
      eval(['schnitzcells(i).' fluor2_backAlt     ' = schnitzcells(i).' fluor2_backAlt_all '( FluorIndex );'])

      eval(['schnitzcells(i).' fluor3_sum         ' = schnitzcells(i).' fluor3_sum_all '( FluorIndex );'])
      eval(['schnitzcells(i).' fluor3_mean        ' = schnitzcells(i).' fluor3_mean_all '( FluorIndex );'])
      eval(['schnitzcells(i).' fluor3_back        ' = schnitzcells(i).' fluor3_back_all '( FluorIndex );'])


  end

   eval(['existmyfield=isfield(schnitzcells, ''' fluor5_mean_all ''');']) % check if field exists (=1)  % take Y5 instead of Y4 for minimalmode
  if existmyfield==1
      if p.minimalMode~=1
        eval(['schnitzcells(i).' fluor4_sum       ' = schnitzcells(i).' fluor4_sum_all '( FluorIndex );'])
        eval(['schnitzcells(i).' fluor4_mean      ' = schnitzcells(i).' fluor4_mean_all '( FluorIndex );'])
      end

    eval(['schnitzcells(i).' fluor5_sum       ' = schnitzcells(i).' fluor5_sum_all '( FluorIndex );'])
    eval(['schnitzcells(i).' fluor5_mean      ' = schnitzcells(i).' fluor5_mean_all '( FluorIndex );'])
    eval(['schnitzcells(i).' fluor5_stdev     ' = schnitzcells(i).' fluor5_stdev_all '( FluorIndex );'])
    eval(['schnitzcells(i).' fluor5_median    ' = schnitzcells(i).' fluor5_median_all '( FluorIndex );'])
  end

  eval(['existmyfield=isfield(schnitzcells, ''' fluor6_mean_all ''');']) % check if field exists (=1)
  if existmyfield==1
    eval(['schnitzcells(i).' fluor6_mean      ' = schnitzcells(i).' fluor6_mean_all '( FluorIndex );'])
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Determine on which fluor fields the next operations will take place
%--------------------------------------------------------------------------
% e.g. fluor2Fields= {'Y' 'Y2' 'Y3'}
if p.minimalMode~=1  
    eval(['fluor2Fields = {''' upper(p.fluorcolor) ''' ''' upper(p.fluorcolor) '2'' ''' upper(p.fluorcolor) '3''};']);
else
    eval(['fluor2Fields = {''' upper(p.fluorcolor) '''};']);
end
eval(['existmyfield=isfield(schnitzcells, ''' fluor5_mean_all ''');']) % check if field exists (=1) %Y5 instead of Y4 for minimalMode
if existmyfield==1
    if p.minimalMode~=1 
        eval(['fluor2Fields(end+1) = {''' upper(p.fluorcolor) '4''};'])
    end
  eval(['fluor2Fields(end+1) = {''' upper(p.fluorcolor) '5''};'])
end
eval(['existmyfield=isfield(schnitzcells, ''' fluor6_mean_all ''');']) % check if field exists (=1)
if existmyfield==1
  eval(['fluor2Fields(end+1) = {''' upper(p.fluorcolor) '6''};'])
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Determine on which fluor fields the next operations will take place
%--------------------------------------------------------------------------
%fluor2Fields = {'R' 'R2' 'R3'};
%if isfield(schnitzcells, 'R4_mean_all')
%  fluor2Fields(end+1) = {'R4'};
%  fluor2Fields(end+1) = {'R5'};
%end
%if isfield(schnitzcells, 'R6_mean_all')
%  fluor2Fields(end+1) = {'R6'};
%end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZCELLS TO ADD FINAL STUFF
%--------------------------------------------------------------------------
% loop over schnitzcells
for i = 1:length(schnitzcells)
  % if exist, get parent, daughters and sister
  clear schnitzP schnitzE schnitzD;
  if (schnitzcells(i).P ~= 0), schnitzP = schnitzcells(schnitzcells(i).P); end
  if (schnitzcells(i).D ~= 0), schnitzD = schnitzcells(schnitzcells(i).D); end
  if (schnitzcells(i).E ~= 0), schnitzE = schnitzcells(schnitzcells(i).E); end

  %------------------------------------------------------------------------
  % Calc average mean fluorescence over cell cycle
  %------------------------------------------------------------------------
  for num = 1:length(fluor2Fields)
    fluor2Field = char(fluor2Fields(num));
    schnitzcells(i).(['av_' fluor2Field '_mean']) = NaN;
    if length(schnitzcells(i).([fluor2Field '_mean'])) > 0
      schnitzcells(i).(['av_' fluor2Field '_mean']) = mean(schnitzcells(i).([fluor2Field '_mean']));
    end
  end
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % Calc the fitted (interpolated) fluor_mean values for each frame
  %------------------------------------------------------------------------
  for num = 1:length(fluor2Fields)
    fluor2Field = char(fluor2Fields(num));

    schnitzcells(i).(['fitted_' fluor2Field '_mean']) = NaN * [1:length(schnitzcells(i).time)];

    % make fluor_time & temp_fluor_mean  (e.g. R_time, temp_r_mean)
    eval([ fluor_time '= [];'])
    eval([ temp_fluor_mean '= [];'])

    if (exist('schnitzP'))
      eval([fluor_time '= [schnitzP.' fluor_time '];'])
      eval([temp_fluor_mean '= [schnitzP.([fluor2Field ''_mean''])];'])
    end

    eval([ fluor_time '= [' fluor_time ' schnitzcells(i).' fluor_time '];'])
    eval([temp_fluor_mean '= [' temp_fluor_mean ' schnitzcells(i).([fluor2Field ''_mean''])];'])

    if (exist('schnitzE') & exist('schnitzD'))
      eval(['data_length         = min( length(schnitzD.' fluor_time '), length(schnitzE.' fluor_time ') );'])
      eval([fluor_time '= [' fluor_time ' schnitzD.' fluor_time '(1:data_length)];'])
      eval([temp_fluor_mean '= [' temp_fluor_mean ' (0.5*(schnitzE.([fluor2Field ''_mean''])(1:data_length) + schnitzD.([fluor2Field ''_mean''])(1:data_length)))];'])
    end

    % interpolate data
    eval(['existlength=length(' fluor_time ');']) % get length (=1)
    if existlength>1
        % sometimes interpolating causes an error (NW2015-09)
        try
          eval(['schnitzcells(i).([''fitted_'' fluor2Field ''_mean'']) = interp1( ' fluor_time ', ' temp_fluor_mean ',  schnitzcells(i).time);'])
        catch
            warning('Could not interpolate fluor field.'); % -MW
            eval(['schnitzcells(i).([''fitted_'' fluor2Field ''_mean'']) = NaN;'])
        end
    end

    % problem with out of range values: are still NaN 
    % replace values by closest known values
    eval([ temp_fitted_fluor_mean '= schnitzcells(i).([''fitted_'' fluor2Field ''_mean'']);'])
    numberLocations = isnan(schnitzcells(i).(['fitted_' fluor2Field '_mean']))-1;
    FirstNumberLocation = find(numberLocations,1,'first');
    LastNumberLocation = find(numberLocations,1,'last');
    if ~isempty(FirstNumberLocation)
      for pos=[1:FirstNumberLocation-1]
        schnitzcells(i).(['fitted_' fluor2Field '_mean'])(pos) = schnitzcells(i).(['fitted_' fluor2Field '_mean'])(FirstNumberLocation);
      end
    end
    if ~isempty(LastNumberLocation)
      for pos=[LastNumberLocation+1,length(schnitzcells(i).(['fitted_' fluor2Field '_mean']))]
        schnitzcells(i).(['fitted_' fluor2Field '_mean'])(pos) = schnitzcells(i).(['fitted_' fluor2Field '_mean'])(LastNumberLocation);
      end
    end

    % Calc the average of fitted fluor_mean values
    schnitzcells(i).(['av_fitted_' fluor2Field '_mean']) = mean(schnitzcells(i).(['fitted_' fluor2Field '_mean']));
  end
  %------------------------------------------------------------------------


  %------------------------------------------------------------------------
%   % Calc the average of fitted R_mean values
%   %------------------------------------------------------------------------
%   for num = 1:length(fluorFields)
%     fluorField = char(fluorFields(num));
% 
%     % problem with NaN values -> replace by closest known values
%     temp_fitted_R_mean = schnitzcells(i).(['fitted_' fluorField '_mean']);
%     numberLocations = isnan(temp_fitted_R_mean)-1;
%     FirstNumberLocation = find(numberLocations,1,'first');
%     LastNumberLocation = find(numberLocations,1,'last');
%     if ~isempty(FirstNumberLocation)
%       for pos=[1:FirstNumberLocation-1]
%         temp_fitted_R_mean(pos) = temp_fitted_R_mean(FirstNumberLocation);
%       end
%     end
%     if ~isempty(LastNumberLocation)
%       for pos=[LastNumberLocation+1,length(temp_fitted_R_mean)]
%         temp_fitted_R_mean(pos) = temp_fitted_R_mean(LastNumberLocation);
%       end
%     end
%     schnitzcells(i).(['av_fitted_' fluorField '_mean']) = mean(temp_fitted_R_mean);
%   end
%   %------------------------------------------------------------------------
end 
%--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % REMOVE FIELDS I WILL PROBABLY NOT USE
% %--------------------------------------------------------------------------
% rmfields = {'autoFluor_all' 'R_sum_all' 'R_mean_all' 'R_back_all' 'R2_sum_all' 'R2_mean_all' ...
%             'R2_back_all' 'R2_backAlt_all' 'R3_sum_all' 'R3_mean_all' 'R3_back_all'};
% schnitzcells = rmfield(schnitzcells, rmfields);
% if isfield(schnitzcells, 'R4_mean_all')
%   rmfields = {'R4_sum_all' 'R4_mean_all' 'R5_sum_all' 'R5_mean_all' 'R5_stdev_all' 'R5_median_all'};
%   schnitzcells = rmfield(schnitzcells, rmfields);
% end
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Display warning if no fluorescence pictures were found
%--------------------------------------------------------------------------
if (fluorcounter==0)
    disp(['Warning! No fluorescence pictures for fluorcolor ''' p.fluorcolor ''' found. ' ...
          'Maybe color in ''p'' wrong associated. Change with p.fluor[x]=...']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Save extacted data
%--------------------------------------------------------------------------
save(p.schnitzName, 'schnitzcells');
disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize performs new normalization
function regNorm = normalize( reg, ...            % gfp image in camera units
                               back, ...           % background in camera units
                               expt, ...           % exposure time in ms
                               binning, ...        % camera bin size (2)
                               micronsPerPixel, ... % micron per pixel (0.04065)
                               autoFluor2)           % autofluorescence of cells
if (expt=='empty')%DE, 11/Nov/2013: was 'emptr', changed to 'empty'.
  disp('Warning! no exposure time specified for fluor image. Taking 1000 arbitrarily.');
  expt = 1000; 
end
regNorm=double(reg)-double(back);
regNorm=regNorm/(expt * binning^2 * micronsPerPixel^2);
regNorm=regNorm-autoFluor2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalizeElowitz performs normalization as Elowitz did
function regNorm = normalizeElowitz( reg, ...     % gfp image in camera units
                                      back, ...    % background in camera units
                                      expt, ...    % exposure time in ms
                                      binning, ... % camera bin size (2)
                                      autoFluor2)    % autofluorescence of cells
regNorm=double(reg)-double(back);
if (expt=='empty')%DE, 11/Nov/2013: was 'emptr', changed to 'empty'.
  disp('Warning! no exposure time specified for fluor image');
  expt = 1; 
end
regNorm=regNorm/expt; % in camera_xfp_1sec_units;
regNorm=regNorm/(binning^2);
regNorm=regNorm-autoFluor2; % in units of molecules?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

