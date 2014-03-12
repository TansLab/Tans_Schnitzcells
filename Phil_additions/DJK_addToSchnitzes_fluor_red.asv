% DJK_addToSchnitzes_fluor_red loads schnitzcells, calculates the fluor for
% schnitzes and saves them to schnitzcells again.
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
%

function [p,schnitzcells] = DJK_addToSchnitzes_fluor_red(p,varargin) 

DEBUG = 1;

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_addToSchnitzes_fluor_red';

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
% Override any schnitzcells parameters/defaults given optional fields/values
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
% If explicit DJK_saveDir2 is not given, define it
if ~existfield(p,'DJK_saveDir2')
  p.DJK_saveDir2 = [p.analysisDir 'fluor2' filesep 'RfromFluor2_' filesep];
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
% Load Schnitzcells
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
trackRange = sort(unique([schnitzcells.frames])) - 1;

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
  for age = 1:length(s.frames)
    framenum = s.frames(age) - 1; % hack - schnitzedit is 1-based, correct here

    % look up that frame number's index within trackRange
    lincellnumIndex = find(trackRange==framenum);
    cellnum  = s.cellno(age);
    lincellnum{lincellnumIndex}(cellnum) = schnitznum;
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER FRAMES IN TRACKRANGE
%--------------------------------------------------------------------------
for num = 1:length(trackRange)
  fr = trackRange(num);

  % empty previous data
  clear Lc rreg rreg4 rreg5 rreg5Norm;

  % load segmented image for this frameNum
  name= [p.segmentationDir,p.movieName,'seg',str3(fr)];
  load([name]); % including variables: LNsub, Lc, phsub, timestamp, rect, yback, yreg, yshift, expty, gainy, rback, rreg, rshift, exptr, gainr    

  % in case onScreen, Lc_block will be used to save where blocks are
  Lc_block = zeros(size(Lc));
  
  %------------------------------------------------------------------------
  % CHECK WHETHER FLUOR EXIST FOR THIS FRAME
  %------------------------------------------------------------------------
  if ~exist('rreg')
    disp([' * Frame ' str3(fr) ' : no fluor2 image']);
    % skip frame if no fluor
  else
    disp([' * Frame ' str3(fr) ' : processing as nr ' str3(num) ' in range']);

    % load fluor for this frameNum
    name= [p.tracksDir, p.movieName, 'Fluor2_', str3(fr)];
    disp(['   -> loading Fluor2 of ', str3(fr)]);
    load([name]); % including variables: exptr, gainr, phaseCropSize, phaseFullSize, rect, rectCrop, rback, rback2, rback2Alt, rbackAlt, rbinning, rreg, rreg2

    %----------------------------------------------------------------------
    % NORMALIZE FLUOR IMAGES
    %----------------------------------------------------------------------
    % rreg    : gedeelte van fluor image, vergroot naar gelang binning
    % rreg2   : flatfield & shading corrected rreg
    % rreg3   : shift corrected rreg2
    % rreg4   : convolved version of rreg2
    % rreg5   : shift corrected rreg4
    % rback   : background of rreg, traditional
    % rback2  : background of rreg2, traditional
    % rback3  : background of rreg3, traditional

    % normalize fluor images: correct for background, exposure time, binning, micronsPerPixel & autoFluor2
    rregNorm  = normalizeElowitz(rreg,rback,exptr,rbinning,p.autoFluor2);
    rreg2Norm = normalize(rreg2,rback2,exptr,rbinning,p.micronsPerPixel,p.autoFluor2);
    rreg3Norm = normalize(rreg3,rback3,exptr,rbinning,p.micronsPerPixel,p.autoFluor2);

    % if deconvolved image exist do the same for them
    if exist('rreg4')
      rreg4Norm = normalize(rreg4,rback2,exptr,rbinning,p.micronsPerPixel,p.autoFluor2);
      rreg5Norm = normalize(rreg5,rback3,exptr,rbinning,p.micronsPerPixel,p.autoFluor2);
    end
    %----------------------------------------------------------------------
  end
  %------------------------------------------------------------------------

  
  %----------------------------------------------------------------------
  % LOOP OVER EACH SCHNITZ THAT EXISTS DURING THIS FRAME, AND UPDATE IT
  %----------------------------------------------------------------------
  schnitzesForFrame = lincellnum{num}; % [schnitznum for cellno1, schnitznum for cellno2, etc]
  nonZeroSchnitzes = schnitzesForFrame(schnitzesForFrame~=0); % with correction you sometimes end up with unexisting schnitzes (0)
  
  for i = nonZeroSchnitzes
    % figure out index within this schnitz' age-based arrays
    age     = find((schnitzcells(i).frames-1) == fr);
    cellno  = schnitzcells(i).cellno(age);
    
    if isempty(age)
      error(['lincellnum says schnitz num ' i 'exists in frame ' fr ', but that frame can''t be found in the schnitz' frames array']);
    end

    schnitzcells(i).autoFluor2_all(age)     = p.autoFluor2;

    %--------------------------------------------------------------------
    % First fill fluor data with NaN
    %--------------------------------------------------------------------
    schnitzcells(i).R_sum_all(age)        = NaN;
    schnitzcells(i).R_mean_all(age)       = NaN;
    schnitzcells(i).R_back_all(age)       = NaN;

    schnitzcells(i).R2_sum_all(age)       = NaN;
    schnitzcells(i).R2_mean_all(age)      = NaN;
    schnitzcells(i).R2_back_all(age)      = NaN;
    schnitzcells(i).R2_backAlt_all(age)   = NaN;

    schnitzcells(i).R3_sum_all(age)       = NaN;
    schnitzcells(i).R3_mean_all(age)      = NaN;
    schnitzcells(i).R3_back_all(age)      = NaN;

    schnitzcells(i).R4_sum_all(age)       = NaN;
    schnitzcells(i).R4_mean_all(age)      = NaN;

    schnitzcells(i).R5_sum_all(age)       = NaN;
    schnitzcells(i).R5_mean_all(age)      = NaN;
    schnitzcells(i).R5_stdev_all(age)     = NaN;
    schnitzcells(i).R5_median_all(age)    = NaN;
    %--------------------------------------------------------------------
    

    %--------------------------------------------------------------------
    % Add fluor data:
    % R_sum_all, SR, R_mean_all, medR, backR, Rall_sum, Rall_stdev, Rall_mean, 
    % Rall_median, Rall_back, Rall_backAlt
    %--------------------------------------------------------------------    
    if exist('rreg')
      
      % get cell number in segmented image for current schnitz & frame
      cellnum = schnitzcells(i).cellno(age);

      % [loc_x, loc_y] are pixels in Lc image where this cell is located
      loc = find(Lc == cellnum); 

      % calc Traditional Fluor data
      schnitzcells(i).R_sum_all(age)        = sum(rregNorm(loc));
      schnitzcells(i).R_mean_all(age)       = schnitzcells(i).R_sum_all(age)/schnitzcells(i).areaPixels(age);
      schnitzcells(i).R_back_all(age)       = rback;

      % calc Flatfield/Shading corrected Fluor data
      schnitzcells(i).R2_sum_all(age)       = sum(rreg2Norm(loc));
      schnitzcells(i).R2_mean_all(age)      = schnitzcells(i).R2_sum_all(age)/schnitzcells(i).areaPixels(age);
      schnitzcells(i).R2_back_all(age)      = rback2;
      schnitzcells(i).R2_backAlt_all(age)   = rback2Alt;

      % calc Flatfield/Shading/Shift corrected Fluor data
      schnitzcells(i).R3_sum_all(age)       = sum(rreg3Norm(loc));
      schnitzcells(i).R3_mean_all(age)      = schnitzcells(i).R3_sum_all(age)/schnitzcells(i).areaPixels(age);
      schnitzcells(i).R3_back_all(age)      = rback3;

      if exist('rreg4')
        % calc Flatfield/Shading/Deconvolved corrected Fluor data
        schnitzcells(i).R4_sum_all(age)     = sum(rreg4Norm(loc));
        schnitzcells(i).R4_mean_all(age)    = schnitzcells(i).R4_sum_all(age)/schnitzcells(i).areaPixels(age);

        schnitzcells(i).R5_sum_all(age)     = sum(rreg5Norm(loc));
        schnitzcells(i).R5_mean_all(age)    = schnitzcells(i).R5_sum_all(age)/schnitzcells(i).areaPixels(age);
        schnitzcells(i).R5_stdev_all(age)   = std(rreg5Norm(loc));
        schnitzcells(i).R5_median_all(age)  = median(rreg5Norm(loc));
      end
    end

    
    %----------------------------------------------------------------------
    % Add block fluor of rreg5Norm
    %----------------------------------------------------------------------
    schnitzcells(i).R6_mean_all(age)      = NaN;
    schnitzcells(i).R6_offset_all(age)    = NaN;

    if exist('rreg5Norm')
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
        offset_measurements(offset_nr).R6_mean = mean( rreg5Norm( offset_measurements(offset_nr).block ) );
%         disp([str3(offset_nr) ' - ' num2str(offset_values(offset_nr)) ' : ' num2str(offset_measurements(offset_nr).R6_mean)]);
      end
      [trash, best_offset_nr] = max( [offset_measurements(:).R6_mean] );

%       disp(['schnitz ' str3(i) ' frame ' str3(fr) ' best_offset_nr ' str3(best_offset_nr) ' best_offset ' num2str(offset_values(offset_nr)) ' value ' num2str(offset_measurements(offset_nr).R6_mean)]);

      schnitzcells(i).R6_mean_all(age)  = offset_measurements(best_offset_nr).R6_mean;
      schnitzcells(i).R6_offset_all(age)    = offset_values(best_offset_nr)* p.micronsPerPixel;

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
  if p.onScreen ~= 2 & exist('rreg5Norm')
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

    % scale yimage
    fluor = round( DJK_scaleRange(rreg5Norm, [0 15], [0 number_of_colors-1]) );
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
    fluor_block( find(Lc_block>0) ) = 1;
    fluor_block = DJK_imagePlace(fluor_block, Lc_bar, 25, 55);
    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------
    % Show figure
    %--------------------------------------------------------------------    
    if p.onScreen == 1
      figureName = [p.movieDate ' ' p.movieName ' fluor2 block fr ' str3(fr)];
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
    if (upper(saveFigInput(1))=='R')
      figureFileName = [p.movieName '-r6-' str3(fr)];
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
  FluorIndex = find(~isnan(schnitzcells(i).R_frames_all));

  schnitzcells(i).autoFluor2       =  schnitzcells(i).autoFluor2_all( FluorIndex );

  schnitzcells(i).R_sum           = schnitzcells(i).R_sum_all( FluorIndex );
  schnitzcells(i).R_mean          = schnitzcells(i).R_mean_all( FluorIndex );
  schnitzcells(i).R_back          = schnitzcells(i).R_back_all( FluorIndex );
  
  schnitzcells(i).R2_sum          = schnitzcells(i).R2_sum_all( FluorIndex );
  schnitzcells(i).R2_mean         = schnitzcells(i).R2_mean_all( FluorIndex );
  schnitzcells(i).R2_back         = schnitzcells(i).R2_back_all( FluorIndex );
  schnitzcells(i).R2_backAlt      = schnitzcells(i).R2_backAlt_all( FluorIndex );

  schnitzcells(i).R3_sum          = schnitzcells(i).R3_sum_all( FluorIndex );
  schnitzcells(i).R3_mean         = schnitzcells(i).R3_mean_all( FluorIndex );
  schnitzcells(i).R3_back         = schnitzcells(i).R3_back_all( FluorIndex );

  if isfield(schnitzcells, 'R4_mean_all')
    schnitzcells(i).R4_sum        = schnitzcells(i).R4_sum_all( FluorIndex );
    schnitzcells(i).R4_mean       = schnitzcells(i).R4_mean_all( FluorIndex );

    schnitzcells(i).R5_sum        = schnitzcells(i).R5_sum_all( FluorIndex );
    schnitzcells(i).R5_mean       = schnitzcells(i).R5_mean_all( FluorIndex );
    schnitzcells(i).R5_stdev      = schnitzcells(i).R5_stdev_all( FluorIndex );
    schnitzcells(i).R5_median     = schnitzcells(i).R5_median_all( FluorIndex );
  end

  if isfield(schnitzcells, 'R6_mean_all')
    schnitzcells(i).R6_mean       = schnitzcells(i).R6_mean_all( FluorIndex );
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Determine on which fluor fields the next operations will take place
%--------------------------------------------------------------------------
fluor2Fields = {'R' 'R2' 'R3'};
if isfield(schnitzcells, 'R4_mean_all')
  fluor2Fields(end+1) = {'R4'};
  fluor2Fields(end+1) = {'R5'};
end
if isfield(schnitzcells, 'R6_mean_all')
  fluor2Fields(end+1) = {'R6'};
end
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
  % Calc the fitted (interpolated) R_mean values for each frame
  %------------------------------------------------------------------------
  for num = 1:length(fluor2Fields)
    fluor2Field = char(fluor2Fields(num));

    schnitzcells(i).(['fitted_' fluor2Field '_mean']) = NaN * [1:length(schnitzcells(i).time)];

    % make R_time & temp_R_mean
    R_time = []; temp_R_mean = [];

    if (exist('schnitzP'))
      R_time = [schnitzP.R_time];
      temp_R_mean = [schnitzP.([fluor2Field '_mean'])];
    end

    R_time = [R_time schnitzcells(i).R_time];
    temp_R_mean = [temp_R_mean schnitzcells(i).([fluor2Field '_mean'])];

    if (exist('schnitzE') & exist('schnitzD'))
      data_length         = min( length(schnitzD.R_time), length(schnitzE.R_time) );
      R_time = [R_time schnitzD.R_time(1:data_length)];
      temp_R_mean = [temp_R_mean (0.5*(schnitzE.([fluor2Field '_mean'])(1:data_length) + schnitzD.([fluor2Field '_mean'])(1:data_length)))];
    end
    
    % interpolate data
    if length(R_time)>1
      schnitzcells(i).(['fitted_' fluor2Field '_mean']) = interp1( R_time, temp_R_mean, schnitzcells(i).time);
    end
    
    % problem with out of range values: are still NaN 
    % replace values by closest known values
    temp_fitted_R_mean = schnitzcells(i).(['fitted_' fluor2Field '_mean']);
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

    % Calc the average of fitted R_mean values
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
% Save extacted data
%--------------------------------------------------------------------------
save(p.schnitzName, 'schnitzcells');
disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize performs new normalization
function rregNorm = normalize( rreg, ...            % gfp image in camera units
                               rback, ...           % background in camera units
                               exptr, ...           % exposure time in ms
                               rbinning, ...        % camera bin size (2)
                               micronsPerPixel, ... % micron per pixel (0.04065)
                               autoFluor2)           % autofluorescence of cells
if (exptr=='emptr')
  disp('Warning! no exposure time specified for fluor image');
  exptr = 1000; 
end
rregNorm=double(rreg)-double(rback);
rregNorm=rregNorm/(exptr * rbinning^2 * micronsPerPixel^2);
rregNorm=rregNorm-autoFluor2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalizeElowitz performs normalization as Elowitz did
function rregNorm = normalizeElowitz( rreg, ...     % gfp image in camera units
                                      rback, ...    % background in camera units
                                      exptr, ...    % exposure time in ms
                                      rbinning, ... % camera bin size (2)
                                      autoFluor2)    % autofluorescence of cells
rregNorm=double(rreg)-double(rback);
if (exptr=='emptr')
  disp('Warning! no exposure time specified for fluor image');
  exptr = 1; 
end
rregNorm=rregNorm/exptr; % in camera_xfp_1sec_units;
rregNorm=rregNorm/(rbinning^2);
rregNorm=rregNorm-autoFluor2; % in units of molecules?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

