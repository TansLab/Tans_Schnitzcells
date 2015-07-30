function p = DJK_makeMovie (p, mode, wrNum, varargin);
% DJK_makeMovie makes movie of cells. Gets data from segmentation files and 
% schnitzcells and can save different kind of images of them. After which a
% movie is made. Full tracking is not required, but 'singleCell' tracking is.
%
% REQUIRED ARGUMENTS:
%
% 'mode' : 'seg'      Color of cells is grey
%          'tree'     Color of cells in subsequent frames are based on tree, such that each schnitz has its own color.
%          'division' Same as tree, but only colors cells just before and after they divide.
%          'phase'    Phasecontrast image. Requires Image Toolbox
%          'xxxxxx'   Color of cells depends on schnitzcells parameter xxxxxx
%
% 'wrNum' : 'none'    No numbers are written
%           'schnitz' Writes schnitz number in cell when cell first occurs
%           'schAll'  Writes schnitz number in every cell every frame
%           'cellno'  Writes cell number in every cell every frame
%           'length'  Draws fit of length in every cell every frame
%
% OPTIONAL ARGUMENTS:
%
% 'outline' = 1       Draws white outline of cells (slow process, default=0)
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
% 'colorNormalize'    In case of fluor mode, the normalization values in
%                     between which the values will be plotted DEFAULT: [0, 12.1]
% 'colorbar' = 0      In case of fluor mode, colorbar can be shown. number indicates number of boxes
%                     (default:10 for fluor, otherwise 0 = not shown)
%
% 'movieOnly' = 1     Assumes images have already been made, and only need
%                     to be sticked together in a movie (default=0)
%
% 'problemCells'      array with problematic cells (+ frame nr), so that
%                     only schnitz nrs will be written in these
%
% 'colorParameter'    In case of fluor mode, datafield from schnitzcells
%                     can be specified here. DEFAULT: 'fitted_Y_mean'
%
% 'schnitzcells'      Use the provided as schnitzcells or lineage data
%                     (instead of loading it from file)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 3;
functionName = 'DJK_makeMovie';

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
case_color = true;
case_fluor = false;
case_tree = false;
case_division = false;
switch upper(mode),
  case 'SEG',
    filnameBase = ['Seg'];
  case 'TREE',
    case_tree = true;
    filnameBase = ['Tree'];
  case 'DIVISION',
    case_division = true;
    filnameBase = ['Division'];
  case 'PHASE',
    case_color = false;
    filnameBase = ['Phase'];
  otherwise % USED TO BE FLUOR MODE
    % Assume that mode is a datafield in schnitzcells
    case_fluor = true;
    p.colorParameter = mode;
    filnameBase = [mode '_'];
end

case_writeNum = true;
case_writeSchnitz = true;
case_writeSchnitzAll = true;
case_drawLength = false;
switch upper(wrNum),
  case 'NONE',
    case_writeNum = false;
  case 'SCHNITZ',
    case_writeSchnitzAll = false;
    filnameBase = [filnameBase 'Num'];
  case 'SCHALL',
    case_writeSchnitzAll = true;
    filnameBase = [filnameBase 'NumAll'];
  case 'CELLNO',
    case_writeSchnitz = false;
    filnameBase = [filnameBase 'Cellno'];
  case 'LENGTH',
    case_writeNum = false;
    case_drawLength = true;
    filnameBase = [filnameBase 'Length'];
  otherwise
    disp(['ERROR: Did not recoginize wrNum input : ' wrNum]);
    return;
end

% Get names of segmentation files in segmentation directory
outprefix = [p.movieName 'seg'];
D = dir([p.segmentationDir, outprefix, '*.mat']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;

% If problemCells mode, numbers are only written for cells/frames in provided problemCells
if ~existfield(p,'problemCells')
    problemCellsOnly = 0;
else
    problemCellsOnly = 1;
end


% If explicit manualRange is not given, take all segmentation files
% Note that the paramter manualRange is used and edited throughout this 
% file to indicate which frames should be analyzed.
% Edit MW 2014/06/24 (see also below)
if ~existfield(p,'manualRange')
    segNameStrings = char(S);
    p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))'; % MW 2014/6/24
    disp('WARNING: p.manualRange was set to contain all segmented frames automatically.')
end

% Determine range which cells are analyzed - edit MW 2014/06/24
if problemCellsOnly % problemcells, select frames with problem cells, also before & after
    errorFramesRange = [];
    errorFrames = unique(p.problemCells(:,2))';
    for fr = errorFrames
      idx = find(p.manualRange == fr);
      if (idx-1>=1)                     , errorFramesRange = [errorFramesRange p.manualRange(idx-1)]; end
      if (idx)                          , errorFramesRange = [errorFramesRange p.manualRange(idx)]; end
      if (idx+1<=length(p.manualRange)) , errorFramesRange = [errorFramesRange p.manualRange(idx+1)]; end
    end
    p.manualRange = unique( errorFramesRange );
end
    
%     fr_min = p.manualRange(1); fr_max = p.manualRange(end); p.manualRange = [];
%     errorFrames = unique(p.problemCells(:,2))';
%     for fr = errorFrames
%       p.manualRange = [p.manualRange fr-1 fr fr+1];
%     end
%     p.manualRange = unique(  p.manualRange(p.manualRange>fr_min-1 & p.manualRange<fr_max+1)  );



disp(['Analyzing ' num2str(length(p.manualRange)) ' frames ', num2str(p.manualRange(1)), ' to ', num2str(p.manualRange(end))]);

% Check whether there are frames without a corrected segmentation
for frameNum = p.manualRange
  clear Lc;
  load([p.segmentationDir, p.movieName, 'seg', str3(frameNum)]);
  if ~(exist('Lc') == 1) 
    disp(['Not corrected segmentation in seg file : ', p.movieName, 'seg', str3(frameNum)]);
  end
end

% If outline, a white outline is added to cells
if ~existfield(p,'outline')
  p.outline = false;
end

% If onScreen, nothing is saved to disc
if ~existfield(p,'onScreen')
  p.onScreen = false;
end

% If stabilize, will try to stabilize image
if ~existfield(p,'stabilize')
  p.stabilize = true;
end

% % If fluor mode, this datafield is used from schnitzcells
% if ~existfield(p,'colorParameter')
%   p.colorParameter = 'fitted_Y_mean';
% end

% If fluor mode, this normalization settings are used
if ~existfield(p,'colorNormalize')
  p.colorNormalize = [0, 12.1];
end

% If colorbar, a colorbar is added on the left top
if ~existfield(p,'colorbar')
  p.colorbar = 0;
  if case_fluor, p.colorbar = 10; end
end

% Set lineage file to load.
if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

% If schnitzcells / lineage file is provided, copy to lineage and remove from p
if existfield(p,'schnitzcells')
  lineage.schnitzcells = p.schnitzcells;
  p = rmfield(p, 'schnitzcells');
end


% if (case_fluor)
%   filnameBase = [p.colorParameter '_' filnameBase];
% end
if (p.outline)
    filnameBase = [filnameBase 'Outline'];
end
if (p.stabilize)
  filnameBase = [filnameBase 'Stabile'];
end
if problemCellsOnly
  filnameBase = [filnameBase 'Problem'];
end

% If explicit DJK_saveDir is not given, define it
if ~existfield(p,'DJK_saveDir')
  datePrefix = [datestr(now,'yyyy-mm-dd') '_' datestr(now,'HHMM') '_'];
  p.DJK_saveDir = [p.analysisDir 'movies' filesep datePrefix filnameBase filesep];
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

% Add folder to filnameBase
filnameBase = [p.DJK_saveDir filnameBase];

% If fluor mode, this datafield is used from schnitzcells
if existfield(p,'movieOnly')
  if p.movieOnly
    % IMMEDIATELY MAKE MOVIES OUT OF THE IMAGES
    DJK_makeMovieFromPNG(p, filnameBase, 'manualRange', p.manualRange);
  end
  return;
end

%disp(filnameBase);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PREPARING TRACKING / SCHNITZCELLS FILE, SLOOKUP AND COLORMAP
%--------------------------------------------------------------------------
if ~exist('lineage') % if schnitzcells provided, no need to load
  %IN CASE OF FLUOR MODE OR LENGTH DRAWING, LOAD SCHNITZCELLS
  if (case_fluor | case_drawLength)
    % Loading an existing schnitz file which contains the image-derived fields
    if exist(p.schnitzName)~=2
      error(['Could not read ' p.schnitzName ' , which is required for fluor mode']);
    else
      lineage = load(p.schnitzName);
      disp(['Load from ''' p.schnitzName ''' completed...']);
    end

    if ( case_fluor & ~isfield(lineage.schnitzcells(1), mode) ) % used to be p.colorParameter
      error(['Seems like ' mode ' is not a field in schnitzcells']); % used to be p.colorParameter
    end

  %IN CASE OF NO FLUOR MODE, LOAD LINEAGE
  else
    % Load lineage file.
    if ~(exist(p.lineageName)==2)
        error(['Could not read lineage file ''' p.lineageName '''.']);
    end
    disp(['Loading lineage from file ''', p.lineageName, '''.']);
    lineage = load(p.lineageName);  
  end
end

% slookup gives the schnitznr for each cellno in each frame.
slookup = makeslookup(lineage.schnitzcells); % note this is local sub func

% Note: in color_map 1 to 251, corresponds to 0 to 250 in unint16 image
% Note: in color_map 1 to 251, corresponds to 1 to 251 in double image
% Make a colormap, will use double image, so image(x,y) = 1, links to color_map(1)
% (will not be used for phase)
if ~case_fluor % number_of_colors is the total number of schnitzes
  number_of_colors = max2(slookup); 
  mymap = DJK_hsv(number_of_colors); 
  [s,I] = sort(rand(number_of_colors,1)); % get sequence of random integers in range
  mymap(2:end+1,:) = mymap(I,:); % randomly reorder mymap color entries, and reserve 1 for background
else % use 250 gradient steps
  number_of_colors = 250;
  mymap = colorGray(number_of_colors); %mymap = jet(number_of_colors); %mymap = ColorSpiral(number_of_colors); %
  mymap(2:end+1,:) = mymap(:,:); % reserve 1 for background
end
% explicitly set some special colors in the colormap
color_background  = 1;
color_text        = number_of_colors+2;
color_outline     = number_of_colors+3;
color_cell        = number_of_colors+4;
mymap(color_background,:) = [0   0   0  ];
mymap(color_text,:)       = [1   1   1  ];
mymap(color_outline,:)    = [0.9 0.9 0.9];
mymap(color_cell,:)       = [0.5 0.5 0.5];

% if colorbar create image that will be used in each frame
if p.colorbar
  % create temporary image containing colorbar
  im_colorbar = zeros(phaseFullSize); 
  % draw colorbar
%   for i=[0:10] % 1 block for p.colorbar values
%     im_colorbar((10*i)+1:(10*(i+1)) , 1:6) = round( DJK_scaleRange(10-i, [0 10], [2 number_of_colors+1]) ); % OLD: M+5+10*(10-i);
  for i=[1:p.colorbar] % 1 block for p.colorbar values
    im_colorbar( 10*i-9:10*i , 1:6 ) = round( DJK_scaleRange(p.colorbar+1-i, [1 p.colorbar], [2 number_of_colors+1]) ); % OLD: M+5+10*(10-i);
  end
  % draw min and max
  min_im = DJK_getImageOfNr( p.colorNormalize(1) );
  max_im = DJK_getImageOfNr( p.colorNormalize(2) );
  min_im(find(min_im>0)) = color_text;
  max_im(find(max_im>0)) = color_text;
  im_colorbar = round( DJK_imagePlace(im_colorbar, min_im, 8+size(min_im,2)/2, p.colorbar*10 - 5) );
  im_colorbar = round( DJK_imagePlace(im_colorbar, max_im, 8+size(max_im,2)/2, 5) );
end
%--------------------------------------------------------------------------


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

  % get the unique nrs in segmentation (uniqueCellNos previously u) 
  % MW 2014/12
  uniqueCellNos = unique(Lc); 
  % get 0 out of the list of unique cellno
  uniqueCellNos = setdiff(uniqueCellNos,0);

  % get the segmented image (subim) 
  subim = double(Lc);
  % subim_color will be the colored image
  subim_color = zeros(size(subim)); % 0 will be background color
  % If number has to be written keep a separate image to put them into
  if (case_writeNum)
    subim_text = zeros(size(subim));
  end

  % If we want a phase image, also get the phase
  if (~case_color)
    phim = double(phsub)/65535;
    phim = imadjust(phim); % adjust the contrast
  end

  %------------------------------------------------------------------------
  % LOOPING OVER CELL NRS
  %------------------------------------------------------------------------
  for uind = 1:length(uniqueCellNos),
    % cellno uind in framenum fr is 'just born', 'growing' or 'dividing'
    temp_schnitz = slookup(fr,uniqueCellNos(uind)); % MW 2014/06/11 removal N+1 bug
    if temp_schnitz==0
      disp(['In frame ' str3(fr) ' cellno ' str3(uind) ' is not linked to schnitz']);
      continue;
    end
    temp_schnitz_P = lineage.schnitzcells(temp_schnitz).P;
    temp_schnitz_frames = lineage.schnitzcells(temp_schnitz).frame_nrs;
    temp_schnitz_cellno = lineage.schnitzcells(temp_schnitz).cellno;
    frIndex = find(temp_schnitz_frames==(fr)); % MW 2014/06/11 removal N+1 bug

    % fl are the locations in image where cellno uind is located
    fl = find(subim==uniqueCellNos(uind));

    %----------------------------------------------------------------------
    % DO THE COLORING
    %----------------------------------------------------------------------
    if (case_color)
      
      % mode=fluor, color cells according to paramater
      if (case_fluor)
        % parameter can exist 1: never
        %                     2: less often then frames
        %                     3: same as frames (can be NaN, than use closest value)
        value = NaN; % in case 1
        if length(lineage.schnitzcells(temp_schnitz).(mode)) > 0
          if ( length(temp_schnitz_frames) > length(lineage.schnitzcells(temp_schnitz).(mode)) )
            % in case 2
            value = lineage.schnitzcells(temp_schnitz).(mode)(1);
          else
            % in case 3
            value = lineage.schnitzcells(temp_schnitz).(mode)(frIndex);
            % in case 3, but NaN
            if isnan(value) 
              not_nan = find(~isnan( lineage.schnitzcells(temp_schnitz).(mode) ));
              if length(not_nan) > 0
                % find closest value
                [min_difference, array_position_in_not_nan] = min(abs(not_nan - frIndex));
                frIndex_closest = not_nan(array_position_in_not_nan);
                value = lineage.schnitzcells(temp_schnitz).(mode)(frIndex_closest);
              end
            end
          end
        end

        % scale
        subim_color(fl) = round(  DJK_scaleRange(value, [p.colorNormalize(1) p.colorNormalize(2)], [2 number_of_colors+1]) );
        % OLD: norm_min = double( p.colorNormalize(1) ); norm_max = double( p.colorNormalize(2) ); 
        % OLD: value = uint16( 0.5 + 99.99*(value-norm_min)/(norm_max-norm_min) ); % normalization
        % OLD: subim_color(fl) = (M+4)+value; % color values are after schnitzcolors
       
      % mode=tree, color cells according to schnitzNr
      elseif (case_tree)
        % color according to schnitzNr+1 (+1 cause else schnitz 1 is black)
        subim_color(fl) = slookup(fr,uniqueCellNos(uind))+1; % MW 2014/06/11 removal N+1 bug
        
      % mode=division, color cells according to schnitzNr, but only at division
      elseif (case_division)
        % just born cells, give color of parent schnitz
        if (temp_schnitz_frames(1)==fr & fr~=1) % MW 2014/06/11 removal N+1 bug
          % in case no parent, give no color
          if temp_schnitz_P>0
            subim_color(fl) = temp_schnitz_P+1;
          else
            subim_color(fl) = color_cell; % light grey
          end

        % cell about to divide, give own color
        elseif (temp_schnitz_frames(end)==fr & fr~=p.manualRange(end)) % MW 2014/06/11 removal N+1 bug
          subim_color(fl) = temp_schnitz+1;
        else
          subim_color(fl) = color_cell; % light grey
        end

      % mode=seg, cells are drawn grey
      else
        % cells are drawn grey
        subim_color(fl) = color_cell; % light grey
      end
    end
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    % ADD OUTLINE
    %----------------------------------------------------------------------
    if (p.outline) 
      test = zeros(size(subim));
      test(fl)=1;
      B = bwboundaries(test, 4);
      B1 = B{1};
      for i = 1:length(B1)
        subim_color(B1(i,1),B1(i,2)) = color_outline;
      end
    end
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    % ADD FITTED LINE
    %----------------------------------------------------------------------
    if (case_drawLength) %& uind<7  & lineage.schnitzcells(temp_schnitz).angle>40 & lineage.schnitzcells(temp_schnitz).angle<50
      % get rotated pixels, so that cell lies straight
      [y,x] = find(subim==uniqueCellNos(uind)); % note: returns (row, column), which will be used as (y,x)
      phi = lineage.schnitzcells(temp_schnitz).rp_angle(frIndex)*pi/180; % convert orientation to radians, and use opposite to rotate back
      x_rot = x*cos(phi) - y*sin(phi); % mathematical rotation
      y_rot = x*sin(phi) + y*cos(phi); % mathematical rotation

      % Get 3rd degree polynomial to rotated pixels
      fitCoef3 = lineage.schnitzcells(temp_schnitz).fitCoef3(frIndex,:);

      % determine for each x in rotated image, what y the fitted line has
      func_3rd = @(x) x.^3 .* fitCoef3(1) + x.^2 .* fitCoef3(2) + x .* fitCoef3(3) + fitCoef3(4);
      
% %       % if bended
% %       if 0
%         x_line_rot = min(x_rot)-2:0.1:max(x_rot)+2;
%         y_line_rot = func_3rd(x_line_rot);
% 
%         % Check which of the rotated pixels are closer than 0.5 away from line
%         clear distances;
%         for i = 1:length(x_rot)
%           distances(i,:) = (x_line_rot - x_rot(i)).^2 + (y_line_rot - y_rot(i)).^2;
%         end
% 
%         idx = find( min(distances,[],2) < 0.5^2);
%         x_rot_intersect = x_rot(idx);
%         y_rot_intersect = y_rot(idx);
% %       else
% %         x_line_rot = min(x_rot)-1:0.5:max(x_rot)+1;
% %         
% %         y_line_rot = func_3rd(x_line_rot]);
% % 
% %         % Check which y_rot is closer than 0.5 away from line
% %         idx_pixel_intersects = find( y_rot>= y_line_rot-0.5 & y_rot<= y_line_rot+0.5 );
% %         x_rot_intersect = x_rot(idx_pixel_intersects);
% %         y_rot_intersect = y_rot(idx_pixel_intersects);
% %       end

    % determine for set of x_rot's, what y_rot the fitted line has
    y_line_rot_left = func_3rd(x_rot-0.5);
    y_line_rot_right = func_3rd(x_rot+0.5);
    for i = 1:length(y_line_rot_left)
      if y_line_rot_left(i) < y_line_rot_right(i)
        temp = y_line_rot_left(i);
        y_line_rot_left(i) = y_line_rot_right(i);
        y_line_rot_right(i) = temp;
      end
    end
    
    % Check which y_rot is closer than 0.5 away from line
%     idx_pixel_intersects = find( y_rot>= y_line_rot-0.5 & y_rot<= y_line_rot+0.5 );
    idx_pixel_intersects = find( y_rot >= y_line_rot_right-0.5 & y_rot <= y_line_rot_left+0.5 );
    x_rot_intersect = x_rot(idx_pixel_intersects);
    y_rot_intersect = y_rot(idx_pixel_intersects);




      
      % Rotate intersect back
      x_intersect = x_rot_intersect*cos(-phi) - y_rot_intersect*sin(-phi); % mathematical rotation
      y_intersect = x_rot_intersect*sin(-phi) + y_rot_intersect*cos(-phi); % mathematical rotation

      rc_intersect = sub2ind(size(subim_color), round(y_intersect), round(x_intersect)) ;
      subim_color( rc_intersect ) = color_outline;
      
%       % New method
%       temp_schnitz_fitCoef3 = lineage.schnitzcells(temp_schnitz).fitCoef3(frIndex,:);
%       temp_schnitz_angle = lineage.schnitzcells(temp_schnitz).angle(frIndex);
%       temp_schnitz_fitXs = lineage.schnitzcells(temp_schnitz).fitXs(frIndex,:);
%       phi = -temp_schnitz_angle*pi/180; % opposite of what is used at fitting
%       fit_x_rot = [temp_schnitz_fitXs(1)-2:0.1:temp_schnitz_fitXs(2)+2]; %thin=:1: thick=:0.0001:
%       fit_y_rot = fit_x_rot .* fit_x_rot .* fit_x_rot .* temp_schnitz_fitCoef3(1) + fit_x_rot .* fit_x_rot .* temp_schnitz_fitCoef3(2) + fit_x_rot .* temp_schnitz_fitCoef3(3) + temp_schnitz_fitCoef3(4);
%       fit_x = round( fit_x_rot*cos(phi) - fit_y_rot*sin(phi) );
%       fit_y = round( fit_x_rot*sin(phi) + fit_y_rot*cos(phi) );
% 
%       fit_rc = sub2ind(size(subim_color), fit_y, fit_x);
%       subim_color( intersect(fit_rc', fl) ) = color_outline;
    end
    %----------------------------------------------------------------------

    
    %----------------------------------------------------------------------
    % WRITE NUMBERS (SCHNITZ NUMBERS OR CELL NUMBERS)
    %----------------------------------------------------------------------
    if (case_writeNum)
      ref_fr = find(temp_schnitz_frames==fr); % MW 2014/06/11 removal N+1 bug

      % if we are in problem mode check whether correct frame/cell, else continue 
      if problemCellsOnly
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
      if (temp_schnitz_frames(1)==fr | case_writeSchnitzAll | ~case_writeSchnitz) % MW 2014/06/11 removal N+1 bug
        subim_text = DJK_imagePlace(subim_text, nr_im, nr_x, nr_y);
      end
    end
    %----------------------------------------------------------------------
    
  end
  %--------------------- END LOOPING OVER CELL NRS  -----------------------

  
  %------------------------------------------------------------------------
  % FINISH IMAGE
  %------------------------------------------------------------------------
  if (case_color) % non-phase image
    if (case_writeNum) % put schnitz/cellno numbers in image
      subim_color(find(subim_text>0)) = color_text;
    end
    
    if (p.onScreen)
      figure('color','k');
      image(ind2rgb(subim_color,mymap));
      axis off;
      axis image;
    else
      % enlarge to original size
      im_color = DJK_enlargeImage(subim_color, phaseFullSize, rect, p.stabilize);
      % put colorbar in left-top
      if p.colorbar 
        im_color(find(im_colorbar>0)) = im_colorbar(find(im_colorbar>0));
      end
      imwrite(ind2rgb(im_color,mymap),[filnameBase str3(fr) '.png'],'png');
      disp(['Writing image ' str3(fr) ': ', [filnameBase str3(fr) '.png']]);
    end

  else % phase image
    if (case_writeNum)
      phim(find(subim_text>0)) = max2(phim);
    end

    % Add outline or length line if they were made
    phim(find(subim_color>0)) = max2(phim);
          
    if (p.onScreen)
      imagesc(phim);
      colormap(gray);
      axis off;
      axis image;
      set(gca,'units','pixels','position',[5 5 size(phim,2)-1 size(phim,1)-1],'visible','off');
    else
      % enlarge to original size
      phim = DJK_enlargeImage(phim, phaseFullSize, rect, p.stabilize);
      imwrite(phim,[filnameBase str3(fr) '.png'],'png');
      disp(['Writing image ' str3(fr) ': ', [filnameBase str3(fr) '.png']]);
    end
  end
end
%--------------------------------------------------------------------------

%keyboard;

%--------------------------------------------------------------------------
% MAKE MOVIES OUT OF THE IMAGES
%--------------------------------------------------------------------------
if (~p.onScreen)
%   DJK_makeMovieFromPNG(p, filnameBase, 'manualRange', p.manualRange); %DJK 090531 disabled, cause never used
end
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makesloopup makes array: schnitz = slookup( frame, cellno) % MW 2014/06/11 removal N+1 bug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slookup = makeslookup(myschnitzcells);
% Creates lookup table that associates frame and cell number to schnitznr

for schnitzIdx = 1:length(myschnitzcells),
  for frameIdx = 1:length(myschnitzcells(schnitzIdx).frame_nrs),
      
    if ~isempty(myschnitzcells(schnitzIdx).frame_nrs),
      if myschnitzcells(schnitzIdx).frame_nrs
          
        slookup(myschnitzcells(schnitzIdx).frame_nrs(frameIdx),myschnitzcells(schnitzIdx).cellno(frameIdx))=schnitzIdx;
        
      end;
    end;
    
  end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
