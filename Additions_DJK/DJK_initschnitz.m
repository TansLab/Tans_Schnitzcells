function p = DJK_initschnitz (movieName, movieDate, movieKind, varargin);
% function P = initschnitz (movieName, movieDate, movieKind, varargin);
%
%   Same as INITSCHNITZ but then allows e.coli.AMOLF movie kind. 
%   
%   -------
%   If no further specs are given will immediately set p.micronsPerPixel = 0.04065
%   (microsocope1, CoolSnap camera)
%   As of 2014-11 a different camera (Hamamatsu) is used on microscope1. For this camera
%   p.micronsPerPixel=0.0434. Call function with
%   (...'camera','hamamatsu'...) to adjust automatically p.micronsPerPixel.
%   ** Once all old experiments are analyzed this micr/px can be set as
%   default **
%   ------
%
%   ----
%   If the microscope2 (new setup) is used initialize with optional
%   argument ('micromanager',1). This will also set the micronsPerPixel to
%   the defaulft 0.04312 [and later enable handling of different data
%   structure from micromanger (vs metamorph)]
%   no further camera specifications are needed.
%   ----
%
%   INITSCHNITZ allows users to get variables and directories set up 
%   properly before performing cell tracking and analysis operations.  
%   It lets users define & control a number of variables (aka parameters) 
%   used during cell tracking and analysis, including a variety of file 
%   names and paths, as well as constants (e.g. min/max cell sizes, kernel 
%   sizes, etc.) that control various phases of the cell segmentation, 
%   tracking and analysis.  It also creates the analysis directories if 
%   necessary.  
%   
%   INITSCHNITZ(movieName,movieDate,movieKind) will initialize and return the 
%   parameter structure P by setting all of the cell analysis parameters to 
%   their default values defined for the given movieKind, and by using default 
%   values for original image and cell analysis paths.
%   
%   INITSCHNITZ(movieName,movieDate,movieKind,'Field1',Value1,'Field2',Value2,...)
%   initializes and returns the parameter structure P by first setting all of 
%   the cell analysis parameters to their default values defined for the given 
%   movieKind, and then overriding specific defaults by setting 
%   P.Field1 = Value1, P.Field2 = Value2, etc.  Thus any/all schnitzcells 
%   parameter values can be defined in the function call via these optional 
%   field/value pairs.  (This is in the style of setting MATLAB properties 
%   using optional property/value pairs.)  
%   
%   The function returns a struct (1x1 struct array) referred to as the 
%   schnitzcells parameter structure that contains fields and values 
%   defined by (or resulting from) the list of properties & values 
%   provided in the optional arguments.
%   
% Required Inputs:
%   
%     movieName           Movie name prefix string used as a portion the 
%                         filenames for the movie in question, e.g. 'movie2-10'
%                         This is also the prefix of all image file names, e.g. 
%                         'movie2-10-p-026.tif'.
%                         
%     movieDate           Date in string format, YYYY-MM-DD, e.g. '2004-04-29'
%                         
%     movieKind           "Kind" of movie, e.g. "bacillus", "e.coli" that 
%                         is used to recall certain analysis parameter values 
%                         that can vary depending on the kind of the movie.
%   
% Directory Structure:
%   
%                              rootDir
%                                 |
%                ------------------------------------------------------------
%                |                                    |                     |
%           2004-04-29                           2004-11-09                ...
%                |                                    |
%        ---------------------------         --------------------------    
%        |         |        |      |         |        |        |      |
%      movie1    movie2   movie3  ...      movie1   movie2   movie3  ...
%        |         |        |                |        |        | 
%       ...        |       ...               |       ...      ...
%                  |                         |
%        ---------------------          ---------------------     
%        |         |         |          |         |         |         
%     images  segmentation tracks    images  segmentation tracks 
% 
%                                
%                                
% Optional Schnitzcells Parameters/Properties Include:
%     
%     rootDir         The top-level directory containing results of cell  
%                     segmentation/tracking analysis for one or more movies.
%                     This directory contains (or INITSCHNITZ will create) 
%                     the YYYY-MM-DD/movieName directory and its 
%                     sub-directories (images,segmentation,tracks).  Under 
%                     the Windows OS this directory path should include the 
%                     drive name, e.g.: 'g:\all-cell-analyses'.  The 
%                     default is simply 'g:\' for Windows, '/' for Unix.
%                     
%     imageDir        You can explicitly specify the directory that contains 
%                     the original images that you want to analyze.  This is 
%                     by default an image directory found under the movie 
%                     directory, but can be overridden by manually specifying 
%                     this variable.  Under the Windows OS this directory path 
%                     should include the drive name, e.g. 'g:\images'
%                     The default is simply 'g:\' for Windows, '/' for Unix.
%                     
%     segRange        A vector of integers specifying the frame numbers 
%                     you intend to segment.  By default this will be all 
%                     frames in the given movie.  Because movies are 
%                     numbered from zero, the values in this range should 
%                     correspond to the numbers of the frames, so 
%                     'segRange',[0:9] will segment the first 10 frames 
%                     numbered 000,001,...,009.
%
%     Fluor1          The fluorescence colour you used. The name of the     % NW (11/12/2)
%                     fluorescence images decides the actual input that is
%                     is necessary. So if e.g. a (red) fluorescence image 
%                     is named pos2-r-001 the input is 'r', if the name 
%                     is pos2-rfp-001, then the input is 'rfp'.
%
%                     ******* Do NEVER use an image with 'none' in its name!!!!
%                     (which should also not be very likely...) *******
%
%     Fluor2, Fluor3  Option to include more fluorescence colors. If none
%                     are used, skip field or assign with 'none' (which
%                     will be done automatically if field was skipped)
%
%     micronsPerPixel Preferably set automatically via micromanager/camera.
%                     Manual setting overwrites the automatic defaults
%
%     The three parameters below will adjust some settings according to 
%     your setup.
%
%     setup           This allows you to let SchnitzCells know which setup
%                     was used for measuring. Settings will be loaded
%                     accordingly. Also defaults for softwarePackage and
%                     camera are adjusted accordingly.
%
%     softwarePackage This allows you to select a software package, i.e.
%                     'metamorph' or 'micromanager' in our case. If
%                     setup=='setup1', 'metamorph' is chosen per default,
%                     if setup=='setup2', micromanager' is chosen by
%                     default.
%                     Note there is also a metadata parameter called 
%                     'Software', which is loaded from the image
%                     information (but this information can only be
%                     extracted with knowledge of which softwarePackage was
%                     used - since metadata is stored differently depending
%                     on the package).
%
%     camera          Important for the micronsPerPixel parameter.
%                     Defualt is 'hamamatsu' when setup=='setup1', and
%                     'unkown' when setup=='setup2'. According to camera
%                     settings, the micronsperpixel will be chosen.
%                     For setup1, available cameras are coolsnap (oldest one) 
%                     and hamamatsu (the new one per ~2014/11/11*). 
%                       *) TODO MW: check this; NW's experiment list has
%                       1st exp. w. new one at 11/11 but might be earlier 
%                       exp. w. new camera
%
%                         
% ----- some old or derived parameters below -----    
%     
%     numphaseslices  now segmoviephase derives this by looking at image dir
%                         
%     

%-------------------------------------------------------------------------------
% Define Constants
%-------------------------------------------------------------------------------

% add any new movie kinds here, and add corresponding logic later in this code
% note movie kinds should be written in lower-case here for string compares
allKnownMovieKinds = {'e.coli','bacillus','e.coli.amolf'}; % DJK 071120

% define every schnitzcells parameter structure field in a bogus structure 
% that we will use at the end simply to impose a field ordering
% --- movie info ---
pOrder.movieName = 0;
pOrder.movieDate = 0;
pOrder.movieKind = 0;
% --- directories ---
pOrder.rootDir = 0;
pOrder.dateDir = 0;
pOrder.movieDir = 0;
pOrder.imageDir = 0;
pOrder.segmentationDir = 0;
pOrder.tracksDir = 0;
pOrder.analysisDir = 0; % DJK 071211
%pOrder.partialDir = 0;
% --- fluorescence colors ---    % NW (11/12/02) (don't know if assignment here is necessary)
pOrder.fluor1=0;
pOrder.fluor2=0;
pOrder.fluor3=0;

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 3;
if (nargin < 3) | ...
   (~isstr(movieName) | ~isstr(movieDate) | ~isstr(movieKind)) | ...
   (mod(nargin,2) == 0) | ...
   (regexp(movieDate, '\d\d\d\d-\d\d-\d\d') ~= 1)
%if ((nargin < 3) | (~isstr(movieName) | ~isstr(movieDate) | ~isstr(movieKind)) | (mod(nargin,2) == 0) | (length(movieDate) ~= 10))
    
    if (regexp(movieDate, '\d\d\d\d-\d\d-\d\d') ~= 1) % MW 2014/12
        disp('*Note that date format wasn''t recognized in movieDate');
    end
    
  errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> initschnitz:',...
      '    Invalid input arguments.',...
      '    Try "help initschnitz".');
  error(errorMessage);
end

% verify that we know what to do with the given movie kind
switch lower(movieKind)
  case allKnownMovieKinds
    % do nothing now; will customize params later
  otherwise
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> initschnitz:',...
        '    Invalid movie kind.',...
        '    Try "help initschnitz", or define new movie kind in initschnitz.');
    error(errorMessage);
end

% quickly iterate over input arguments just to do a little testing
% this may not be necessary because we loop over varargin this way later on
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  % loop over pairs of remaining/optional arguments
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
          'Error using ==> initschnitz:',...
          '    Invalid property ', num2str(varargin{i}), ...
	  ' is not (needs to be) a string.',...
          '    Try "help initschnitz".');
      error(errorMessage);
    end
    % OK, I lied, I'm doing a *little* more than just input testing...
    % if rootDir given in input args, grab it now (ignore others until later)
    fieldName = DJK_schnitzfield(varargin{i});
    switch (fieldName)
      case 'rootDir'
%        p.(fieldName) = varargin{i+1}; this is how general case was coded
        p.rootDir = varargin{i+1};
    end
  end
end


%-------------------------------------------------------------------------------
% Fill the schnitzcells parameter structure given required inputs and defaults
%-------------------------------------------------------------------------------

% move the main inputs into the parameter structure
% (This *could* be overwritten if user again specifies any of 
%  movieName,movieDate,movieKind passed in optional args!)
p.movieName = movieName;
p.movieDate = movieDate(1:10); % MW 2014/12
p.movieKind = movieKind;

% define default analysis root directory
if (~existfield(p,'rootDir'))
  p.rootDir = 'g:\';
else
  % strip any trailing slash because below code assumes it's not there 
  pathlength = length(p.rootDir);
  if (p.rootDir(pathlength) == filesep)
    p.rootDir = p.rootDir(1:end-1);
  end
end

% set up variables that hold date information 
movie_year  = movieDate(1:4);
movie_month = movieDate(6:7);
movie_day   = movieDate(9:10);
datelong = str2num([movie_year movie_month movie_day]);

fs = filesep;

p.dateDir         = [p.rootDir,fs,movieDate,fs];
p.movieDir        = [p.rootDir,fs,movieDate,fs,movieName,fs];
p.imageDir        = [p.rootDir,fs,movieDate,fs,movieName,fs,'images',fs];
p.segmentationDir = [p.rootDir,fs,movieDate,fs,movieName,fs,'segmentation',fs];
p.tracksDir       = [p.rootDir,fs,movieDate,fs,movieName,fs,'data',fs];
p.analysisDir     = [p.rootDir,fs,movieDate,fs,movieName,fs,'analysis',fs]; % DJK 071211
%p.partialDir      = [p.rootDir,fs,movieDate,fs,movieName,fs,'partial',fs];

%-------------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values 
% to the schnitzcells parameter structure
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
          'Error using ==> initschnitz:',...
          '    Invalid property ', num2str(varargin{i}), ...
	  ' is not (needs to be) a string.',...
          '    Try "help initschnitz".');
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
    % also need to append placeholder for any custom fields to pOrder structure
    if (~existfield(pOrder,fieldName))
      pOrder.(fieldName) = 0; % place holder for ordering
    end
  end
end

% MW 2014/11
% if setup is not set manually, assume setup1 
if ~existfield(p,'setup') 
    p.setup='setup1';
    disp('WARNING: no setup specified, will assume you used microscope 1!');
    pOrder.setup = 0; % place holder for ordering
end

% if software package is not set manually, asume Metamorph (= setup1)
% note there is also 'p.Software', which is actually the software version 
% (This is a TODO.)
if ~existfield(p,'softwarePackage') 
    if strcmp(p.setup,'setup1')
        p.softwarePackage='metamorph'; % default software for setup1
    elseif strcmp(p.setup,'setup2') 
        p.softwarePackage='micromanager'; % default software for setup2
    else
        p.softwarePackage='metamorph'; % default software in general
    end
    disp(['WARNING: no software Package (metamorph or micromanager) specified, will assume you used ' p.softwarePackage '!']);
    pOrder.softwarePackage = 0; % place holder for ordering
end

% if camera is not set manually, set it now to hamamatsu (the setup1 new one from 2014/11)
if ~existfield(p,'camera')
    if strcmp(p.setup,'setup1')
        p.camera='coolsnap'; % default camera for setup1 % MW 2014/12 backwards compatib. (= old camera not used any more)
    elseif strcmp(p.setup,'setup2')
        p.camera='hamamatsu2'; % default camera for setup2
    else
        p.camera='unknown'; % default camera in general
    end
    pOrder.camera=0; % place holder for ordering
    disp(['WARNING: no camera defined, assuming ' p.camera ' !']);
end

% If micronsPerPixel has not been set manually, set it now.
% Depends on setup.
if ~existfield(p,'micronsPerPixel')
    switch p.setup  % default values depend on used camera
        case 'setup1'
            if strcmp(p.camera, 'hamamatsu')
                % hamamatsu camera, setup 1 (new camera)
                p.micronsPerPixel = 0.0438;
                    % (installed in 2014-11)
            elseif strcmp(p.camera, 'coolsnap')
                % CoolSnap camera, setup 1 (old camera)
                p.micronsPerPixel = 0.04065;
            end
        case 'setup2'
            % camera hamamatsu, setup 2 (always same camera)
            p.micronsPerPixel = 0.04312;
        otherwise
            error('Couldn''t set micronsPerPixel')
    end
    pOrder.micronsPerPixel = 0;
end 

% If less than 3 fluorescence colours have been set, fill others with
% 'none'                                                                 NW 11/12/02
if ~existfield(p,'fluor1')
    p.fluor1='none';
end
if ~existfield(p,'fluor2')
    p.fluor2='none';
end
if ~existfield(p,'fluor3')
    p.fluor3='none';
end


%-------------------------------------------------------------------------------
% if a fillinator 'cell' optional argument is provided, adjust the movie 
% subdirectory names accordingly
%-------------------------------------------------------------------------------
if existfield(p,'cell')

  % adjust segmentation dir (after checking for ending filesep)
  pathlength = length(p.segmentationDir);
  if (p.segmentationDir(pathlength) == filesep)
    p.segmentationDir = p.segmentationDir(1:end-1);
  end
  p.segmentationDir = [p.segmentationDir '-cell-' str3(p.cell) filesep];

  % adjust tracks dir (after checking for ending filesep)
  pathlength = length(p.tracksDir);
  if (p.tracksDir(pathlength) == filesep)
    p.tracksDir = p.tracksDir(1:end-1);
  end
  p.tracksDir = [p.tracksDir '-cell-' str3(p.cell) filesep];

%   % adjust partial dir (after checking for ending filesep)
%   pathlength = length(p.partialDir);
%   if (p.partialDir(pathlength) == filesep)
%     p.partialDir = p.partialDir(1:end-1);
%   end
%   p.partialDir = [p.partialDir '-cell-' str3(p.cell) filesep];

  % define a default fillinatorFile
  if ~existfield(p,'fillinatorFile')
    p.fillinatorFile = [p.segmentationDir,'cell-',str3(p.cell),'-pixels.mat'];
  end

  % add extra fields to pOrder since these new fillinator fields were added
  pOrder.cell = 0;
  pOrder.fillinatorFile = 0;

end


%-------------------------------------------------------------------------------
% Act on the variables now that they're all defined (make directories & such)
%-------------------------------------------------------------------------------

% check if image directory already exists; use info to warn later on
image_dir_exists_before = false;
if exist(p.imageDir) == 7
  image_dir_exists_before = true;
end

% if directories don't exist, create them  (will issue warnings if they exist)
if exist(p.dateDir)~=7
  [status,msg,id] = mymkdir(p.dateDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.dateDir ' : ' msg]);
  end
end
if exist(p.movieDir)~=7
  [status,msg,id] = mymkdir(p.movieDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.movieDir ' : ' msg]);
  end
end
if exist(p.imageDir)~=7
  [status,msg,id] = mymkdir(p.imageDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.imageDir ' : ' msg]);
  end
end
if exist(p.segmentationDir)~=7
  [status,msg,id] = mymkdir(p.segmentationDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.segmentationDir ' : ' msg]);
  end
end
if exist(p.tracksDir)~=7
  [status,msg,id] = mymkdir(p.tracksDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.tracksDir ' : ' msg]);
  end
end
if exist(p.analysisDir)~=7 % DJK 071211
  [status,msg,id] = mymkdir(p.analysisDir); % DJK 071211
  if status == 0 % DJK 071211
    disp(['Warning: unable to mkdir ' p.analysisDir ' : ' msg]); % DJK 071211
  end % DJK 071211
end % DJK 071211


% if image directory did not exist before we created it, warn user that 
% segmentation won't work until it's populated
if ~image_dir_exists_before & exist(p.imageDir)==7
  disp(['Warning: The image directory ' p.imageDir ' was just created.']);
  disp(['         You will need to move images there before segmentation.']);
end

% make sure every directory field has a trailing filesep
pathlength = length(p.rootDir);
if (p.rootDir(pathlength) ~= filesep)
  p.rootDir = [p.rootDir filesep];
end
pathlength = length(p.dateDir);
if (p.dateDir(pathlength) ~= filesep)
  p.dateDir = [p.dateDir filesep];
end
pathlength = length(p.movieDir);
if (p.movieDir(pathlength) ~= filesep)
  p.movieDir = [p.movieDir filesep];
end
pathlength = length(p.imageDir);
if (p.imageDir(pathlength) ~= filesep)
  p.imageDir = [p.imageDir filesep];
end
pathlength = length(p.segmentationDir);
if (p.segmentationDir(pathlength) ~= filesep)
  p.segmentationDir = [p.segmentationDir filesep];
end
pathlength = length(p.tracksDir);
if (p.tracksDir(pathlength) ~= filesep)
  p.tracksDir = [p.tracksDir filesep];
end
pathlength = length(p.analysisDir); %DJK 071211
if (p.analysisDir(pathlength) ~= filesep) %DJK 071211
  p.analysisDir = [p.tracksDir filesep]; %DJK 071211
end %DJK 071211


% make sure we order the fields properly before returning
p = orderfields(p, pOrder);

return

