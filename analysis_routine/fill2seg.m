
function s = fill2seg(p,varargin)

% FILL2SEG  Convert fillinator cells/track to schnitzcells segmentation files
%
%   FILL2SEG allows users to convert a fillinator output file containing 
%   one cell's pixel locations over a range of frames to a series of 
%   schnitzcells-compliant segmentation files.  This permits the fillinator 
%   results to be reviewed/adjusted with manualcheckseg, and subsequently 
%   fed into the existing fluorescence analysis pipeline (i.e. schnitzabsaugen, 
%   plotschnitzme, etc.). 
%
%   FILL2SEG reads from a file known as the fillinatorFile that contains a 
%   struct array containing frame, xpixels, ypixels for each frame processed. 
%   A segmentation file is produced for each frame in the fillinatorFile.
%
%   FILL2SEG(P,'Field1',Value1,'Field2',Value2,...) does the same thing,
%   but the extra arguments permit users to adjust any parameters controlling
%   the process.  The extra arguments can override any specific parameters
%   provided in P by setting P.Field1 = Value1, P.Field2 = Value2, etc.
%   Thus any/all schnitzcells parameter values can be defined in the
%   function call via these optional field/value pairs.  (This is in the
%   style of setting MATLAB properties using optional property/value pairs.)
%
%   FILL2SEG returns a struct (1x1 struct array) referred to as the
%   schnitzcells parameter structure that contains fields and values
%   contained in P, including unchanged/original parameters plus any of those
%   added or overridden in the list of properties & values provided in the
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control FILL2SEG
%
%   fillinatorFile filename of fillinator pixel data file generated (default
%                  filename is [p.segmentationDir 'cell-001-pixels.mat'])
%
%   regsize        the maximum fluorescent image offset used by quickreg
%
%   offByOne       (HACK) a flag that when set permits dealing with old 
%                  fillinator results that have 1-based frame numbers, 
%                  rather than the default 0-based frame numbers.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) || ...
		(mod(nargin,2) == 0) || ...
		(~isSchnitzParamStruct(p))
	errorMessage = sprintf ('%s\n%s\n%s\n',...
		'Error using ==> fill2seg:',...
		'    Invalid input arguments.',...
		'    Try "help fill2seg".');
	error(errorMessage);
end

%-------------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
	for i=1:2:(numExtraArgs-1)
		if (~isstr(varargin{i}))
			errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
				'Error using ==> fill2seg:',...
				'    Invalid property ', num2str(varargin{i}), ...
				' is not (needs to be) a string.',...
				'    Try "help fill2seg".');
			error(errorMessage);
		end
		fieldName = schnitzfield(varargin{i});
		p.(fieldName) = varargin{i+1};
	end
end

if ~existfield(p,'regsize')
  % maximum value of any translation between phase and fluorescent images
  p.regsize = 3;
end

% An optional offByOne argument (flag) handles older off-by-one fillinator data
if ~existfield(p,'offByOne')
    p.offByOne = false
end

% Load the fillinator results file 
if ~existfield(p,'fillinatorFile')
  p.fillinatorFile = [p.segmentationDir,'cell-001-pixels.mat'];
end

if ~exist(p.fillinatorFile) == 2
  error(['Could not read fillinator file ''' p.fillinatorFile '''.']);
end

[pathstr,name,ext,versn] = fileparts(p.fillinatorFile);
disp (['Converting fillinator output ' name ' to segmentation files...']);
load(p.fillinatorFile)

% figure out phase image filename convention (so we can trust numphaseslices)
% Get phase image(s)
mnamePhase2   = [p.movieName,'-p-2-*.tif'];
mnamePhaseAll = [p.movieName,'-p-*.tif'];

Dphase2   = dir([p.imageDir, mnamePhase2]);
DphaseAll = dir([p.imageDir, mnamePhaseAll]);

if p.numphaseslices == 1
  D = DphaseAll;
  if isempty(D)
    errorMessage = sprintf('%s%s\n%s%s\n',...
        '    No images found in directory ', p.imageDir, ...
        '    matching pattern ', mnamePhaseAll);
    error(errorMessage);
  end
else
  D = Dphase2;
  if isempty(D)
    errorMessage = sprintf('%s%s\n%s%s\n',...
        '    No images found in directory ', p.imageDir, ...
        '    matching pattern ', mnamePhase2);
    error(errorMessage);
  end
end


allFrames = sort([fillstruc.frame]);

% Loop over all frames in the fillinator output
for frameNum = allFrames 

  frameNumString = str3(frameNum);
  disp(['Creating segmentation output for frame ' frameNumString '...']);

  % First, load the fillinator results for this frame, and double-check them
  currFillstruc = fillstruc(find([fillstruc.frame] == frameNum));
  % JCR HACK: set p.offByOne to use old fillinator data with one-too-high frames
  if p.offByOne
    currFillstruc = fillstruc(find([fillstruc.frame] == frameNum+1));
  end
  if isempty(currFillstruc)
    disp(['Frame ' frameNumString ' has no cell pixels, skipping.']);
    continue
  end
  if isempty(currFillstruc.xpixels)
    disp(['Frame ' frameNumString ' has no cell pixels, skipping.']);
    continue
  end

  % a segmented frame should have at least these fields, so start making them
  savelist=['''phsub'',''LNsub'',''LC'',''rect'',''timestamp'''];

  % Load the phase image(s) into phsub
  if (p.numphaseslices == 1)
    phaseImgName = [p.imageDir p.movieName '-p-' frameNumString '.tif'];
    phaseImg = imread(phaseImgName);
  else
    phaseImgName = [p.imageDir p.movieName '-p-2-' frameNumString '.tif'];
    phaseImg = imread(phaseImgName);
  end
  fullsize = [size(phaseImg,1) size(phaseImg,2)];

  % create phsub given fillinator output pixels
  extra= 75; % <- extra number of pixels on either side of region of interest
  xmin= max(min(currFillstruc.xpixels) - extra, 1);
  xmax= min(max(currFillstruc.xpixels) + extra, size(phaseImg,1));
  ymin= max(min(currFillstruc.ypixels) - extra, 1);
  ymax= min(max(currFillstruc.ypixels) + extra, size(phaseImg,2));
  % In spite of code in segmoviephase suggesting rect = [xmin ymin xmax ymax];
  % rect is truely [rowmin rowmax colmin colmax] = [ymin ymax xmin xmax];
  % Likewise, fillinator's xpixels are row indices, and ypixels are col indices.
  % So while all this is WRONG, at least it's all consistent
  rect = [ xmin ymin xmax ymax ];
  phsub= phaseImg(xmin:xmax, ymin:ymax);

  % create LNsub given fillinator output pixels
  LNsub = zeros(size(phsub,1),size(phsub,2));
  for i=1:length(currFillstruc.xpixels)
    LNsub(currFillstruc.xpixels(i)-xmin,...
          currFillstruc.ypixels(i)-ymin) = 1;
  end
  
  % create LC also to indicate this segmentation is already corrected
  LC = LNsub;

  [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings(phaseImgName,'p');

  % get fluorescence files and registration (when available)
  cname= [p.imageDir,p.movieName,'-c-',frameNumString,'.tif'];
  yname= [p.imageDir,p.movieName,'-y-',frameNumString,'.tif'];
        
  if exist(cname)==2
    disp('found CFP image');
    [creg, cshift, cback]= quicknoreg(LNsub, cname, rect, p.regsize, fullsize);
    [exptcstr, gainc, exptc]= imsettings(cname);
    savelist=[savelist,',''creg'',''cshift'',''exptc'',''gainc'',''cback'''];
  end
  if exist(yname)==2
    disp('found YFP image');
    [yreg, yshift, yback]= quicknoreg(LNsub, yname, rect, p.regsize, fullsize);	
    [exptystr, gainy, expty]= imsettings(yname);
    savelist=[savelist,',''yreg'',''yshift'',''expty'',''gainy'',''yback'''];
  end
          
  segmentationOutFile = [p.segmentationDir p.movieName 'seg' frameNumString];

  if exist([segmentationOutFile,'.mat']),
      disp('appending...');
      savelist=[savelist,',''-append'''];
  end;
  
  % have to save via eval to be flexible for varying fluorescence files
  eval(['save(''',segmentationOutFile,''',',savelist,');']);

end % for frameNum = allFrames

