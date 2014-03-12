function p = DJK_segmoviephaseImproved (p,varargin)
% DJK_SEGMOVIEPHASEIMPROVED Segment movie from phase images.
%
% copied from DJK_SEGMOVIEPHASE, which was copied from SEGMOVIEPHASE with
% the following changes: 
%
% * Only has p.movieKind option 'e.coli.AMOLF' with different paramater
%   values than e.coli and bacillus
%
%   Use of DJK_segphaseImproved() instead of DJK_segphase().
%
%   Removed p.segmentationPhaseSlice and p.prettyPhaseSlice
%   Addition of p.edgeSlices and p.fillingEdgeSlices
%   
%   if paramater 'DJK_saveDir' is provided, segmentation and intermediat L
%   images are also saved to this subdirectory of the segmentation
%   directory. 
%
%   Parameter 'noKinky' is added. If it is true, then the step of breaking
%   up kinky cells is skipped. Default: noKinky is true.
%
%   Parameter 'maxCellWidthConservative' is added. If it is < maxCellWidth
%   than before cutting cells between L5 and L6, first a cut is done with a 
%   smaller maxCellWidth. It is set to 6 for e.coli.AMOLF, maxCellWidth otherwise.
%
%   useMedfilt2forEdge
%   useL1MaskCutting
%
%
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control segmoviephase
% 
%   segRange                range of frame numbers to segment; by default 
%                           all image frames in imageDir will be segmented 
%   
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
      'Error using ==> segmoviephase:',...
      '    Invalid input arguments.',...
      '    Try "help segmoviephase".');
  error(errorMessage);
end
%-------------------------------------------------------------------------------


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
          'Error using ==> segmoviephase:',...
          '    Invalid property ', num2str(varargin{i}), ...
          ' is not (needs to be) a string.',...
          '    Try "help segmoviephase".');
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i}); % DJK 071122
    p.(fieldName) = varargin{i+1};
  end
end

% DJK_saveDir indicates where is provided, image
if isfield(p,'DJK_saveDir')
  % make sure every directory field has a trailing filesep
  pathlength = length(p.DJK_saveDir);
  if (p.DJK_saveDir(pathlength) ~= filesep)
    p.DJK_saveDir = [p.DJK_saveDir filesep];
  end
  % if directory doesn't exist, create it
  if exist(p.DJK_saveDir)~=7
    [status,msg,id] = mymkdir([p.segmentationDir,p.DJK_saveDir]);
    if status == 0
      disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    end
  end
  % Als make sure png directory exists
  png_dir = [p.segmentationDir 'png' filesep]; 
  if exist(png_dir)~=7
    [status,msg,id] = mymkdir([png_dir]);
    if status == 0
      disp(['Warning: unable to mkdir ' png_dir ' : ' msg]);
      return;
    end
  end
else
  p.DJK_saveDir = [''];
end

% figure out (automatically) if we have just 1 phase or multiple phase images
% This affects/selects which phase images to use... at least initially 
% for determining number and range of frames in the movie.
mnamePhase2   = [p.movieName,'-p-2-*.tif'];
mnamePhaseAll = [p.movieName,'-p-*.tif'];
Dphase2   = dir([p.imageDir, mnamePhase2]);
DphaseAll = dir([p.imageDir, mnamePhaseAll]);

if ~isfield(p,'numphaseslices')
  D = Dphase2; % for segRange later
  % Derive the number of phase slices by examining the image directory files.
  % We understand 3 cases:
  %  1) if one   phase slice,  Dphase2 is empty
  %  2) if two   phase slices, DphseAll is 2x size of Dphase2
  %  3) if three phase slices, DphseAll is 3x size of Dphase2
  if isempty(DphaseAll)
    error(['Can''t find any images in directory ' p.imageDir]);
  end
  if isempty(Dphase2)
    D = DphaseAll; % for segRange later
    p.numphaseslices = 1;
    disp('You appear to have 1 phase image per frame: setting p.numphaseslices = 1;');
  elseif length(DphaseAll) == length(Dphase2)*2 
    p.numphaseslices = 2;
    disp('You appear to have 2 phase images per frame: setting p.numphaseslices = 2;');
  elseif length(DphaseAll) == length(Dphase2)*3 
    p.numphaseslices = 3;
    disp('You appear to have 3 phase images per frame: setting p.numphaseslices = 3;');
  elseif length(DphaseAll) == length(Dphase2)*4 %DJK 080803
    p.numphaseslices = 4; %DJK 080803
    disp('You appear to have 4 phase images per frame: setting p.numphaseslices = 4;'); %DJK 080803
  elseif length(DphaseAll) == length(Dphase2)*5 %DJK 080214
    p.numphaseslices = 5; %DJK 080214
    disp('You appear to have 5 phase images per frame: setting p.numphaseslices = 5;'); %DJK 080214
  else
    error(['Can''t figure out how many phase slices are in your movie! (' num2str(length(DphaseAll)) ',' num2str(length(Dphase2)) ')']);
  end
end

% Define segmentation parameters that vary with Movie Kind here:
switch lower(p.movieKind)
  case 'e.coli.amolf'
    disp('Movie kind: e.coli.AMOLF');
    
    % STEP 0A: FILTER AND EDGE PHASE
    if ~isfield(p,'edgeSlices') % DJK 081228
      p.edgeSlices = [ 1:p.numphaseslices ]; % slices that are used for edge detection
      %'e.coli' p.edgeSlices = [ p.segmentationPhaseSlice ]; %'bacillus' p.edgeSlices = [ p.segmentationPhaseSlice ];
    end
    if ~isfield(p,'useMedfilt2forEdge') % DJK 081228
      p.useMedfilt2forEdge = 0; % whether medfilt2 is applied on images for edge 
      %'e.coli' p.useMedfilt2forEdge = 1; %'bacillus' p.useMedfilt2forEdge = 1;
    end
    if ~isfield(p,'edge_lapofgauss_sigma') % CHECKED -> 2 ok
      p.edge_lapofgauss_sigma = 2; % DJK 081211
      %'e.coli' p.edge_lapofgauss_sigma = 2; %'bacillus' p.edge_lapofgauss_sigma = 3;
    end

    % STEP 0B: FIND MASK_1 (smallest region containing cells)
    if ~isfield(p,'minNumEdgePixels')
      p.minNumEdgePixels = 250; % DJK 071122
      %'e.coli' p.minNumEdgePixels = 250; %'bacillus' p.minNumEdgePixels = 215;  
    end

    % STEP 0D: EXTRACT SMALLER SUBREGIONS

    % STEP 0E: MAKE AN EVEN SMALLER MASK: MASK_2

    % STEP 1: FILL CELLS 
    if ~isfield(p,'L1Fill_only_use_mask_2_bigger') % DJK 090108
      p.L1Fill_only_use_mask_2_bigger = 0
      ; % option to only use mask_2_bigger, and not mask_1 (default: 0, original: 1 & 0)
      %'e.coli' p.L1Fill_only_use_mask_2_bigger = 1; %'bacillus' p.L1Fill_only_use_mask_2_bigger = 1;
    end

    % STEP 2: USING BOTTOM HAT, GET SEEDS TO REMOVE REGION BETWEEN CELLS
    if ~isfield(p,'fillingEdgeSlices') % DJK 081228
      p.fillingEdgeSlices = [ 1:p.numphaseslices ];
      %'e.coli' p.fillingEdgeSlices = [ p.segmentationPhaseSlice ]; %'bacillus' p.fillingEdgeSlices = [ p.segmentationPhaseSlice ];
    end
    if ~isfield(p,'botHatSize') % DJK 081229
      p.botHatSize = 5; % bigger botHatSize gives bigger seeds, sometimes 6 seems to work better
      %'e.coli' p.botHatSize = 5; %'bacillus' p.botHatSize = 5;
    end
    if ~isfield(p,'L2_increase_mask_2_bigger') % DJK 090614
      p.L2_increase_mask_2_bigger = 0;
    end
    
    % STEP 3: CLEAN UP UNSIGHTLY CELLS
    if ~isfield(p,'minCellLengthConservative')
      p.minCellLengthConservative = 20; % DJK 071127: 20 DJK 090531: 30
      %'e.coli' p.minCellLengthConservative = 20; %'bacillus' p.minCellLengthConservative = 12;
    end
    
    % STEP 4: FILL IN EDGE (using seeds from L3)
    
    % STEP 5A: CUTTING BY REMOVING SINGLE PIXELS THAT CONNECT DISPARATE CELLS

    % STEP 5B: CUTTING AT BRANCHPOINTS USING cutcurvbpts
    if ~isfield(p,'maxThreshCut')
      p.maxThreshCut = 0.3;
      %'e.coli' p.maxThreshCut = 0.3; %'bacillus' p.maxThreshCut = 0.3;
    end
    if ~isfield(p,'maxCellWidth')
      p.maxCellWidth = 8; % DJK 071127
      %'e.coli' p.maxCellWidth = 7; %'bacillus' p.maxCellWidth = 20;
    end    
    
    % STEP 5C: REMOVE SMALL CELLS

    % STEP 6A: CUTTING CELLS AT CONCAVE POINTS USING cutcurv
    if ~isfield(p,'waistCut') % DJK 090105
      p.waistCut = 1; % When 0, turns off cutting at waist in 6A & 6B, default: turned on
      %'e.coli' p.waistCut = 1; %'bacillus' p.waistCut = 1;
    end
    if ~isfield(p,'maxThreshCut2')
      p.maxThreshCut2 = 0.2; % smaller maxThreshCut2 the more points are cut 
      %'e.coli' p.maxThreshCut2 = 0.2; %'bacillus' p.maxThreshCut2 = 0.2;
    end
    if ~isfield(p,'maxCellWidthConservative')
      p.maxCellWidthConservative = 8; % two points must be closer than this for cut to be accepted
      %'e.coli' p.maxCellWidthConservative = p.maxCellWidth; %'bacillus' p.maxCellWidthConservative = p.maxCellWidth;
    end 
    if ~isfield(p,'minCellLength')
      p.minCellLength = 30; % cut ignored if it creates `cell' smaller than this
      %'e.coli' p.minCellLength = 30; %'bacillus' p.minCellLength = 20;
    end    
    
    % STEP 6B: CUTTING CELLS AT CONCAVE POINTS USING cutcurv

    % STEP 7: CUT LONG CELLS BY LOOKING AT PHASECONTRAST USING breakcell
    if ~isfield(p,'maxThresh')
      p.maxThresh = 0.25;
      %'e.coli' p.maxThresh = 0.25; %'bacillus' p.maxThresh = 0.25;
    end   
    if ~isfield(p,'minThresh')
      p.minThresh = 0.25;
      %'e.coli' p.minThresh = 0.25; %'bacillus' p.minThresh = 0.25;
    end
    
    % STEP 9A: REMOVE SMALL CELLS
    if ~isfield(p,'removeSmallCells') % DJK 090114
      p.removeSmallCells = 1; % Turns on removing small cells at end, default: turned on, original: on
      %'e.coli' p.removeSmallCells = 1; %'bacillus' p.removeSmallCells = 1;
    end
    
    % STEP 9B: INCREASE CELL WITH 1 PIXEL SO EDGE PIXELS ARE INCLUDED

    % STEP 9C: RESIZE OUTPUT IMAGES TO EXACTLY FIT SEGMENTED CELLS

  otherwise
    errorMessage = sprintf ('%s\n%s\n',...
        'Error using ==> segmoviephase:',...
        '    Invalid movie kind setting internal movie kind parameters.');
    error(errorMessage);
end



% if user didn't specify a range of frames to segment, figure it out 
% given the image names
if ~isfield(p,'segRange')
  [s,I] = sort({D.name}');
  D = D(I);
  numpos= findstr(D(1).name, '.tif')-3;
  imageNameStrings = char(s);
  p.segRange = str2num(imageNameStrings(:,numpos:numpos+2))';
  clear s I numpos imageNameStrings;
end

% Display settings that will be used
disp ('using schnitzcells parameter structure:');
disp (p);
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% Loop over frames in range: will load images and do segmentation
%-------------------------------------------------------------------------------
for i= p.segRange

  % Load images into X
  Dframe = dir([p.imageDir p.movieName '*-p*-' str3(i) '.tif']); %list of all phase contrast images of this frame in directory
  pname = Dframe(1).name; % first filename
  if p.numphaseslices==1
    ph3(:,:,1) = imread([p.imageDir,pname]);
    disp(['reading ',p.imageDir,pname]);
  else
    fstr = findstr(pname,'-p-');
    for islice = 1:p.numphaseslices
      pname(fstr+3) = num2str(islice);
      ph3(:,:,islice) = imread([p.imageDir,pname]);
      disp(['reading ',p.imageDir,pname,' as slice ',num2str(islice)]);
    end
  end

  % Do segmentation on X
  % ph3   : the phase contrast images of this frame in full size
  % phsub : the phase contrast images of this frame in size that contains segmented cells
  % LNsub : the segmented image of this frame in size that contains segmented cells
  % rect  : transformation required to reconstruct full size image from smaller ones
  % phseg : the phase contrast image that was used for edge detection (can be average of multiple phase contrast images)
  [phsub,LNsub,rect,phseg]= DJK_PN_segphaseImproved(ph3, p);

  % DJK 081229 Used to save phsub(:,:,p.prettyPhaseSlice) as phsub
  % From now on will save phseg as phsub
  phsub = phseg;
  
  % Prepare seg data for saving
  if isempty(LNsub) %if no cells segmented, add this so that we have an image
    LNsub = [ 0 ]; 
  end 
  
  [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings([p.imageDir,Dframe(1).name],'p'); % DJK -> will always use timestamp of first phase image
  clear TST*
  
  phaseFullSize = size( ph3(:,:,1) ); % size of full image
  
  savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'''];

  % Prepare fluor data for saving
  yname = [p.imageDir,p.movieName,'-y-',str3(i),'.tif'];
  if exist(yname)==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error
    disp('found Fluor image');
    [yreg, yshift, yback, ybinning] = quicknoreg(LNsub,yname,rect,0,phaseFullSize); % not really important, cause redone during DJK_correctFluorImage
    [exptystr, gainy, expty] = imsettings(yname);
    savelist=[savelist,',''yreg'',''expty'',''gainy'',''yback'',''ybinning'''];
  end

  % Save segmentation file
  if isfield(p,'DJK_saveDir')
    Lname = [p.DJK_saveDir, p.movieName, 'seg', str3(i)];   
  else
    Lname = [p.movieName, 'seg', str3(i)];   
  end
  eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
  disp(['saved file ',p.segmentationDir,Lname]);

  % Save extra seg.png image
  if isfield(p,'DJK_saveDir')
    edgeSlices = num2str(p.edgeSlices);
    edgeSlices(edgeSlices==' ') = [];
    fillingEdgeSlices = num2str(p.fillingEdgeSlices);
    fillingEdgeSlices(fillingEdgeSlices==' ') = [];
    extra_string_filename = '';
    if p.removeSmallCells==0, extra_string_filename = [extra_string_filename '_removeSmallCells0']; end
    png_filename = [png_dir 'seg' str3(i) '_improved_edge' edgeSlices '_fill' fillingEdgeSlices '_pix' str3(p.minNumEdgePixels) '_hat' num2str(p.botHatSize) extra_string_filename '.png']; % '_removeSmallCells0' 
    % DJK_writeSegImage(LNsub, png_filename);
    imwrite(DJK_imshowlabel(LNsub,'phaseImage',phseg), png_filename, 'png','bitdepth',8); % DJK 090529 get seg on top of phase
  end
end
%-------------------------------------------------------------------------------


