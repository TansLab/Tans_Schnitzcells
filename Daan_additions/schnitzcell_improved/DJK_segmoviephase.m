function p = DJK_segmoviephase (p,varargin)
% DJK_SEGMOVIEPHASE   Segment movie from phase images.
% 
%   Same as SEGMOVIEPHASE, but then with p.movieKind option 'e.coli.AMOLF' 
%   that has different paramater values:
%
%   Use of DJK_segphase() instead of segphase().
%
%   Corrected error when working with p.segmentationPhaseSlice = 2 and
%   p.prettyPhaseSlice = 1
%   
%   if paramater 'DJK_saveDir' is provided, segmentation and intermediat L
%   images are also saved to this subdirectory of the segmentation
%   directory. 
%
%   Parameter 'noKinky' is added. If it is true, then the step of breaking
%   up kinky cells is skipped. If noKinky is not added, it is false.
%
%   Parameter 'useBigMask' is added. If it is true, then at 1 point in the 
%   analysis, bigmask is used instead of mask, to avoid having missing cells.
%   It is false for e.coli.AMOLF, false otherwise.
%
%   Parameter 'maxCellWidthConservative' is added. If it is < maxCellWidth
%   than before cutting cells between L5 and L6, first a cut is done with a 
%   smaller maxCellWidth. It is set to 6 for e.coli.AMOLF, maxCellWidth otherwise.
%
%   SEGMOVIEPHASE(P) segments the movie from phase images according to the 
%   parameters and controls described in the schnitzcells parameter structure P,
%   typically generated using INITSCHNITZ.
%   
%   SEGMOVIEPHASE(P,'Field1',Value1,'Field2',Value2,...) segments the movie 
%   from phase images according to the parameters and controls described in the 
%   schnitzcells parameter structure P, but only after P has been updated by 
%   setting P.Field1 = Value1, P.Field2 = Value2, etc.  Thus any schnitzcells 
%   parameters can be updated (overridden) in the function call via these 
%   optional field/value pairs.  (This is in the style of setting MATLAB 
%   properties using optional property/value pairs.)  For a complete list of 
%   schnitzcells parameter structure fields, see INITSCHNITZ.
%   
%   SEGMOVIEPHASE returns a struct array, the updated schnitzcells parameter 
%   structure p that reflects fields updated via any optional field/value 
%   arguments, as well as fields updated to reflect the status of the 
%   segmentation at the time the segmentation exited.
%   
%   For example, P = SEGMOVIEPHASE(P,'segRange',[5:10]) would perform 
%   segmentation on frames 005 thru 010 of the movie described in P, and 
%   return the updated schnitzcells parameter structure P.
%   
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control segmoviephase
% 
%   segRange                range of frame numbers to segment; by default 
%                           all image frames in imageDir will be segmented 
%   
%   prettyPhaseSlice        the phase image number with best focus; with more 
%                           than one phase slice this defaults to 2, otherwise
%                           defaults to 1
%                           
%   segmentationPhaseSlice  the phase image number used when performing cell 
%                           boundary segmentation; defaults to 1
%   
%-------------------------------------------------------------------------------
%

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

% START DJK 071120
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
    % START DJK 081021    
    % Als make sure png directory exists
    png_dir = [p.segmentationDir 'png' filesep]; 
    if exist(png_dir)~=7
      [status,msg,id] = mymkdir([png_dir]);
      if status == 0
        disp(['Warning: unable to mkdir ' png_dir ' : ' msg]);
        return;
      end
    end
    % END DJK 081021
end
% END DJK 071120

% Define segmentation parameters that vary with Movie Kind here:
switch lower(p.movieKind)
   
% START DJK 071120
  case 'e.coli.amolf'
    disp('Movie kind: e.coli.AMOLF');
    if ~isfield(p,'edge_lapofgauss_sigma')
      p.edge_lapofgauss_sigma = 2;
    end
    if ~isfield(p,'minCellArea')
      %p.minCellArea = 300; % DJK 071122
      p.minCellArea = 150; % DJK 081211 -> cause of cell missing sometimes
    end
    if ~isfield(p,'minCellLengthConservative')
      p.minCellLengthConservative = 20; % DJK 071127
    end
    if ~isfield(p,'minCellLength')
      p.minCellLength = 45; % DJK 071122
    end
    if ~isfield(p,'maxCellWidth')
      p.maxCellWidth = 8; % DJK 071127
    end
    if ~isfield(p,'minNumEdgePixels')
      p.minNumEdgePixels = 250; % DJK 071122
    end
    if ~isfield(p,'maxThreshCut')
      p.maxThreshCut = 0.3;
    end
    if ~isfield(p,'maxThreshCut2')
      p.maxThreshCut2 = 0.2;
    end
    if ~isfield(p,'maxThresh')
      p.maxThresh = 0.25;
    end
    if ~isfield(p,'minThresh')
      p.minThresh = 0.25;
    end
    if ~isfield(p,'imNumber1')
      p.imNumber1 = 2;
    end
    if ~isfield(p,'imNumber2')
      p.imNumber2 = 1;
    end
    if ~isfield(p,'radius')
      p.radius = 5;
    end
    if ~isfield(p,'angThresh')
      p.angThresh = 2.7;
    end
    if ~isfield(p,'useBigMask') %DJK 071126
      p.useBigMask = false; %DJK 071126
    end %DJK 071126
    if ~isfield(p,'maxCellWidthConservative') %DJK 071127
      p.maxCellWidthConservative = 6; %DJK 071127
    end %DJK 071127
    if ~isfield(p,'noKinky') %DJK 071127
      p.noKinky = false; %DJK 071127
    end %DJK 071127
% END DJK 071120

  case 'e.coli'
    disp('Movie kind: e.coli');
    if ~isfield(p,'edge_lapofgauss_sigma')
      p.edge_lapofgauss_sigma = 2;
    end
    if ~isfield(p,'minCellArea')
      p.minCellArea = 30;
    end
    if ~isfield(p,'minCellLengthConservative')
      p.minCellLengthConservative = 20;
    end
    if ~isfield(p,'minCellLength')
      p.minCellLength = 30;
    end
    if ~isfield(p,'maxCellWidth')
      p.maxCellWidth = 7;
    end
    if ~isfield(p,'minNumEdgePixels')
%     p.minNumEdgePixels = 300;  % JCR: was dropped to 250 from original 300
     p.minNumEdgePixels = 250;  % JCR: value used for nature methods paper
    end
    if ~isfield(p,'maxThreshCut')
      p.maxThreshCut = 0.3;
    end
    if ~isfield(p,'maxThreshCut2')
      p.maxThreshCut2 = 0.2;
    end
    if ~isfield(p,'maxThresh')
      p.maxThresh = 0.25;
    end
    if ~isfield(p,'minThresh')
      p.minThresh = 0.25;
    end
    if ~isfield(p,'imNumber1')
      p.imNumber1 = 2;
    end
    if ~isfield(p,'imNumber2')
      p.imNumber2 = 1;
    end
    if ~isfield(p,'radius')
      p.radius = 5;
    end
    if ~isfield(p,'angThresh')
      p.angThresh = 2.7;
    end
    if ~isfield(p,'useBigMask') %DJK 071126
      p.useBigMask = false; %DJK 071126
    end %DJK 071126
    if ~isfield(p,'maxCellWidthConservative') %DJK 071127
      p.maxCellWidthConservative = p.maxCellWidth; %DJK 071127
    end %DJK 071127
    if ~isfield(p,'noKinky') %DJK 071127
      p.noKinky = false; %DJK 071127
    end %DJK 071127

  case 'bacillus'
    disp('Movie kind: bacillus');
    if ~isfield(p,'edge')
      p.edge_lapofgauss_sigma = 3;
    end
    if ~isfield(p,'minCellArea')
      p.minCellArea = 30;
    end
    if ~isfield(p,'minCellLengthConservative')
      p.minCellLengthConservative= 12; % JCR observed 2005-02-16
    end
    if ~isfield(p,'minCellLength')
      p.minCellLength = 20;            % JCR observed 2005-02-16
    end
    if ~isfield(p,'maxCellWidth')
      p.maxCellWidth = 20;             % JCR observed 2005-02-16
    end
    if ~isfield(p,'minNumEdgePixels')
      p.minNumEdgePixels = 215;        % ME set to 215 in mail to JCR 2005-08-18
    end
    if ~isfield(p,'maxThreshCut')
      p.maxThreshCut = 0.3;
    end
    if ~isfield(p,'maxThreshCut2')
      p.maxThreshCut2 = 0.2;
    end
    if ~isfield(p,'maxThresh')
      p.maxThresh = 0.25;
    end
    if ~isfield(p,'minThresh')
      p.minThresh = 0.25;
    end
    if ~isfield(p,'imNumber1')
      p.imNumber1 = 2;
    end
    if ~isfield(p,'imNumber2')
      p.imNumber2 = 1;
    end
    if ~isfield(p,'radius')
      p.radius = 5;
    end
    if ~isfield(p,'angThresh')
      p.angThresh = 2.7;
    end
    if ~isfield(p,'useBigMask') %DJK 071126
      p.useBigMask = false; %DJK 071126
    end %DJK 071126
    if ~isfield(p,'maxCellWidthConservative') %DJK 071127
      p.maxCellWidthConservative = p.maxCellWidth; %DJK 071127
    end %DJK 071127
    if ~isfield(p,'noKinky') %DJK 071127
      p.noKinky = false; %DJK 071127
    end %DJK 071127
    
  otherwise
    errorMessage = sprintf ('%s\n%s\n',...
        'Error using ==> segmoviephase:',...
        '    Invalid movie kind setting internal movie kind parameters.');
    error(errorMessage);
end

disp ('using schnitzcells parameter structure:');
disp (p);

% figure out (automatically) if we have just 1 phase or multiple phase images
% This affects/selects which phase images to use... at least initially 
% for determining number and range of frames in the movie.
mnamePhase2   = [p.movieName,'-p-2-*.tif'];
mnamePhaseAll = [p.movieName,'-p-*.tif'];

Dphase2   = dir([p.imageDir, mnamePhase2]);
DphaseAll = dir([p.imageDir, mnamePhaseAll]);


if ~isfield(p,'numphaseslices')
  % Derive the number of phase slices by examining the image directory files.
  % We understand 3 cases:
  %  1) if one   phase slice,  Dphase2 is empty
  %  2) if two   phase slices, DphseAll is 2x size of Dphase2
  %  3) if three phase slices, DphseAll is 3x size of Dphase2
  if isempty(DphaseAll)
    error(['Can''t find any images in directory ' p.imageDir]);
  end
  if isempty(Dphase2)
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

% This code is not robust: problem if there are 2 or 3 phases and you set 
% numphaseslices = 1, D is too big.  Smarter to figure out which D is 
% smaller, and use the smaller one  (BTW, same code's in expand_segimage)
if p.numphaseslices == 1
  D = DphaseAll;
  if length(D) == 0
    errorMessage = sprintf('%s%s\n%s%s\n',...
        '    No images found in directory ', p.imageDir, ...
        '    matching pattern ', mnamePhaseAll);
    error(errorMessage);
  end
else
  D = Dphase2;
  if length(D) == 0
    errorMessage = sprintf('%s%s\n%s%s\n',...
        '    No images found in directory ', p.imageDir, ...
        '    matching pattern ', mnamePhase2);
    error(errorMessage);
  end
end

% Look at phase image file names to figure out how many frames we have
[s,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.tif')-3;

% if user didn't specify a range of frames to segment, figure it out 
% given the image names
if ~isfield(p,'segRange')
  imageNameStrings = char(s);
  p.segRange = str2num(imageNameStrings(:,numpos:numpos+2))';
end

% set some defaults if they don't exist
if exist('outprefix')~=1 
  outprefix = [p.movieName 'seg'];
end
if ~isfield(p,'prettyPhaseSlice');
    if p.numphaseslices>1
        p.prettyPhaseSlice = 2;
    else
        p.prettyPhaseSlice = 1;
    end
end
if ~isfield(p,'segmentationPhaseSlice');
    p.segmentationPhaseSlice = 1;
end
regsize= 0;% max translation (pixels) between phase and fl. images - obsolete.
%%%ST HERE THE PC IMAGES ARE READ
%%%-------------------------------------------------
for i= p.segRange

%    disp(['Daan: p.prettyPhaseSlice = ',int2str(p.prettyPhaseSlice)]);
%    disp(['Daan: p.segmentationPhaseSlice = ', int2str(p.segmentationPhaseSlice)]);
%    disp(['Daan: p.numphaseslices = ', int2str(p.numphaseslices)]);
    
    mynum = str3(i);
    Dframe = dir([p.imageDir p.movieName '*-p*-' str3(i) '.tif']);

    pname = Dframe(1).name;
    if p.numphaseslices==1
        X(:,:,1) = imread([p.imageDir,pname]);
        disp(['reading ',p.imageDir,pname]); % moved - DJK 071127
    else
        fstr = findstr(pname,'-p-');
        for islice = 1:p.numphaseslices
            pname(fstr+3) = num2str(islice);
            X(:,:,islice) = imread([p.imageDir,pname]);
            disp(['reading ',p.imageDir,pname,' as slice ',num2str(islice)]); % extra - DJK 071127
        end
    end

    %[LN,phsub,LNsub,rect,L8,L7,L6,L5,L4]= DJK_segphaseTestAverageImages(X, p); %DJK 081211     
    [LN,phsub,LNsub,rect,L8,L7,L6,L5,L4]= DJK_segphase(X, p); %DJK 071017
    
    % rect
    phaseFullSize = size(LN);
    
    % START DJK 071120
    % Changed because gave error when p.segmentationPhaseSlice = 2 and
    % p.prettyPhaseSlice = 1
    if p.segmentationPhaseSlice~=p.prettyPhaseSlice
      phseg = phsub(:,:,p.segmentationPhaseSlice);
    end
    if size(phsub,3)>1, % nitzan June25: this used to be: phsub = phsub(:,:,2);
      phsub = phsub(:,:,p.prettyPhaseSlice); 
    end;
    savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'''];
    if p.segmentationPhaseSlice~=p.prettyPhaseSlice
      savelist = [savelist,',''phseg'''];
    end
    % END DJK 071120
        
    [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings([p.imageDir,Dframe(1).name],'p'); % DJK -> will use timestamp of first phase image
    clear TST*

    Lname= [outprefix,mynum];   
    cname= [p.imageDir,p.movieName,'-c-',mynum,'.tif'];
    yname= [p.imageDir,p.movieName,'-y-',mynum,'.tif'];
    gname= [p.imageDir,p.movieName,'-g-',mynum,'.tif'];
    rname= [p.imageDir,p.movieName,'-r-',mynum,'.tif'];
    if exist(cname)==2
      disp('found CFP image');
      [creg, cshift, cback, cbinning]= quicknoreg(LNsub,cname,rect,regsize,phaseFullSize);
      [exptcstr, gainc, exptc]= imsettings(cname);
      savelist=[savelist,',''creg'',''cshift'',''exptc'',''gainc'',''cback'',''cbinning'''];
    end
    if exist(yname)==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error DJK 081105
      disp('found YFP image');
      %keyboard
      [yreg, yshift, yback, ybinning]= quicknoreg(LNsub,yname,rect,regsize,phaseFullSize);
      [exptystr, gainy, expty]= imsettings(yname);
      savelist=[savelist,',''yreg'',''yshift'',''expty'',''gainy'',''yback'',''ybinning'''];
      %keyboard
    end
      
    if exist(gname)==2
      disp('found GFP image');
      [greg, gshift, gback, gbinning]= quicknoreg(LNsub,gname,rect,regsize,phaseFullSize);
      [exptgstr, gaing, exptg]= imsettings(gname);
      savelist=[savelist,',''greg'',''gshift'',''exptg'',''gaing'',''gback'',''gbinning'''];
    end
      
    if exist(rname)==2
      disp('found RFP image');
      [rreg, rshift, rback, rbinning]= quicknoreg(LNsub,rname,rect,regsize,phaseFullSize);
      [exptrstr, gainr, exptr]= imsettings(rname);
      savelist=[savelist,',''rreg'',''rshift'',''exptr'',''gainr'',''rback'',''rbinning'''];
    end
      
    % START DJK 071120
    if isfield(p,'DJK_saveDir')
        if isempty(LNsub) %DJK 080219
            LNsub = [ 0 ]; %DJK 080219
        end %DJK 080219
        DJK_writeSegImage(LNsub,[p.segmentationDir p.DJK_saveDir 'seg.png']);
        
        %DJK_writeSegImage(LNsub,[p.segmentationDir,DJK_convertSaveDir(p.DJK_saveDir),'.png']); % DJK 081021
        png_filename = [png_dir 'seg' str3(i) '_slice' num2str(p.segmentationPhaseSlice) '_pixels' str3(p.minNumEdgePixels) '.png']; % DJK 081021
        DJK_writeSegImage(LNsub, png_filename); % DJK 081021
        
        eval(['save(''',p.segmentationDir,p.DJK_saveDir,Lname,''',',savelist,');']);
        disp(['saved file ',p.segmentationDir,p.DJK_saveDir,Lname]);
    else
        eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
        disp(['saved file ',p.segmentationDir,Lname]);
    end
    %eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
    %disp(['saved file ',p.segmentationDir,Lname]);
    % END DJK 071120
    
    clear L*;
end
