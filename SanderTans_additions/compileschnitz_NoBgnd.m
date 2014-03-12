function [p,schnitzcells] = compileschnitz_NoBgnd(p,varargin) 
% COMPILESCHNITZ  given a lineage, suck out fluor info & calc area, length, etc
% 
%   COMPILESCHNITZ allows users to extract fluorescence info from a 
%   cell lineage tree that describes how cells are tracked across frames of 
%   a movie.  The fluorescence info is attached to the tree, and the resulting 
%   annotated tree, called the "schnitz", is written to the schnitzName file.
%   
%   COMPILESCHNITZ(P,'Field1',Value1,'Field2',Value2,...) also performs 
%   cell fluorescence extraction, but the extra arguments permit users 
%   to adjust any parameters describing the movie or parameters controlling 
%   the fluorescence extraction process.  The extra arguments can override 
%   any specific parameters provided in P by setting P.Field1 = Value1, 
%   P.Field2 = Value2, etc.  Thus any/all schnitzcells parameter values can 
%   be defined in the function call via these optional field/value pairs.  
%   (This is in the style of setting MATLAB properties using optional 
%   property/value pairs.)  
%   
%   COMPILESCHNITZ returns a struct (1x1 struct array) referred to as the 
%   schnitzcells parameter structure that contains fields and values 
%   contained in P, including unchanged/original parameters plus any of those 
%   added or overridden in the list of properties & values provided in the 
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters you can adjust to control COMPILESCHNITZ:
% 
%   lineageName    name of un-annotated input lineage name (default filename is
%                  [movieName '_lin.mat'])
%                  
%   schnitzName    name of annotated output schnitz file (default filename is 
%                  [movieName '-Schnitz.mat'])
%                  
%   flatFieldName  name of file with flat field parameters
%                  
%   quickMode      flag, if true, will make it run in quick mode, i.e. it 
%                  will skip extracting the fluorescence from the images 
%                  and only calculate derived fields.  Default value is false.
%   
%   load           flag, if true, will load the saves schnitz file. 
%
%   autoCFL        default = 0
%   yfpOffset      default = 0.
%   
%   crosstalkRatioC2Y    default = 0.0005
%   micronsPerPixel      default = 1/15
%   
%-------------------------------------------------------------------------------
%

%-------------------------------------------------------------------------------
% Variable renaming during rewrite:
%   dir_images  -> p.imageDir
%   schnitzname -> p.schnitzName
%   lineagename -> p.lineageName
%   magic       -> p.magic        (Michael says default = 1, users can override)
%   beadCFP0510 -> p.beadCFP
%   beadYFP0510 -> p.beadYFP
%   Auto_CFL    -> p.autoCFL
%   Auto_YFL    -> p.autoYFL
%   Auto_GFL    -> p.autoGFL
%   Auto_RFL    -> p.autoRFL
%   microns_per_pixel -> micronsPerPixel
%   crosstalk_ratio_C2Y -> crosstalkRatioC2y

% variables needed by old code:
%   calibration
%-------------------------------------------------------------------------------

% old inputs: 
%   catalognumber
% older inputs:
%   framerange_unshifted
%   cellmat
%   If no inputs are given, data is taken from namesdata.mat file.


% creates cell structure from cellmat and images - nitzan's recipe. 
%   P:          parent
%   S:          sister
%   D:          daughter (containing oldest pole)
%   E:          other daughter
%   B:          birth frame
%   N:          life time
%   hist:       lineage within schnitzcells
%   frames:     frames cells existed in
%   cellno:     cell number in image and in cellmat
%   SZ:         cell size (area)
%   ang:        orientation (from regionprops) (given in degrees)
%   wid:        cell width (minoraxislength)
%   len:        cell length (majoraxislength)
%   width_Microns:        cell width (minoraxislength) in Microns
%   lengthMicrons:        cell length (majoraxislength) in Microns
%   volume:     our best estimate of the volume of the cell in micron^3
%   cenx:       x coordinate of centroid in LNsub
%   ceny:       y coordinate of centroid in LNsub
%   cenX:       x coordinate of centroid in full image
%   cenY:       y coordinate of centroid in full image
%   thetas:     angle of orientation axis (given in radians)
%                   defined from old pole (at origin) to new pole
%   datetime:       actual date/time that the (phase) image was taken -- this
%                     is a matlab "datenumber"
%   mins:       number of minutes elapsed since the first frame in the movie.
%   gen:        number of generations since beginning of movie (first cells are gen=0)
%   gens:       "gen" field of length = # of frames
%   phase:      how far is the cell into its cell-cycle (current frame / total frames)
%
%   The following are created for each color of fluorescence (C -> Y -> G -> R etc.)
%   The first three have NaN for frames without fluorescence images in this channel:
%   FC:         total fluorescence intensity (for CFP)
%   MC:         mean fluorescence intensity (for CFP)
%   SC:         standard deviation of fluorescence intensity (for CFP)
%   The rest have values only for frames with fluorescence images in this channel,
%   their length is equal to the number of values which aren't NaN in the above:
%   FCs:        total fluorescence intensity (for CFP)
%   MCs:        mean fluorescence intensity (for CFP)
%   LCs:        fluorescence intensity (CFP) divided by (length*12) (12=average width)
%   FCsframes:  frames which had fluorescence images in this channel (CFP)
%   FCsmins:    times of fluorescence data (CFP)
%   dFCdt:      time derivative of flurescence (dFC/dt)
%   dFCmins:    time of derivatives, = mid-point between the times of the two FCs.


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
   (mod(nargin,2) == 0) | ...
   (~isSchnitzParamStruct(p))
  errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> trackcomplete:',...
      '    Invalid input arguments.',...
      '    Try "help trackcomplete".');
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
          'Error using ==> absaugen:',...
          '    Invalid property ', num2str(varargin{i}), ...
          ' is not (needs to be) a string.',...
          '    Try "help absaugen".');
      error(errorMessage);
    end
    fieldName = schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end

% Set default parameter values if they don't exist yet

if ~existfield(p,'lineageName')
  p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end

if ~existfield(p,'flatFieldName')
  p.flatFieldName = [p.dateDir , 'M22calibrations.mat'];
end%  p.flatFieldName = [p.imageDir,p.movieName,'bead_disk_calibrations.mat'];

if ~existfield(p,'quickMode')
  p.quickMode = 0;
end

if ~existfield(p,'camera_bin_size')
  p.camera_bin_size = 1;
  disp('setting: p.camera_bin_size = 1;');
end
if ~existfield(p,'micronsPerPixel')
  p.micronsPerPixel = 1/15;
end
if ~existfield(p,'autoCFL')
  p.autoCFL = 0;
end
if ~existfield(p,'yfpOffset')
  p.yfpOffset = 0;
end
if ~existfield(p,'autoGFL')
  p.autoGFL = 0;
end
if ~existfield(p,'autoRFL')
  p.autoRFL = 0;
end
if ~existfield(p,'crosstalkRatioC2Y')
  p.crosstalkRatioC2Y = 0; % if images are uncalibrated best not make a mess!
end
if existfield(p,'load') & p.load==1 & exist(p.schnitzName)==2
  disp(['Loading schnitz file ',p.schnitzName,' (p.load=1)']);
  load(p.schnitzName, 'schnitzcells');
  return;
end
  

% Begin the work of updating the schnitz's information

if ~p.quickMode
    % add fields to schnitz lineage by sucking fluorescence info from images
    % requires a lineage, segmentation and images.
    if ~(exist(p.lineageName)==2)
        error(['Could not read lineage file ''' p.lineageName '''.']);
    end
    disp(['Loading lineage from file ''' p.lineageName ...
          ''' and extracting image-based fields...']);
    load(p.lineageName);

    % We need trackRange- we can derive it 
    trackRange = sort(unique([schnitzcells.frames])) - 1; % JCR: hack, 1-based!
                                                          % trackRange 0-based!
    
    % For convenience/efficiency, first create the lincellnum structure
    % initialize lincellnum to have zero'd arrays for each frame 
    lincellnum = {};
    for lincellnumIndex = 1:length(trackRange)
        segdata = load([p.segmentationDir,p.movieName,'seg',...
                        str3(trackRange(lincellnumIndex))]);
        if ~isfield(segdata,'Lc');      % nitzan added 2005June25
            segdata.Lc = segdata.LNsub; % nitzan added 2005June25
        end                             % nitzan added 2005June25
        lincellnum{lincellnumIndex} = zeros([1 max2(segdata.Lc)]);
    end
    % Now, step through each schnitz, store schnitz number in lincellnum
    for schnitznum = 1:length(schnitzcells)
        s = schnitzcells(schnitznum);
        schnitzcells(schnitznum).B = min(s.frames); % just for convenience later
        schnitzcells(schnitznum).N = length(s.frames); % nitzan 2005May16th
        schnitzcells(schnitznum).len = []; % JCR 2005-06-23
        schnitzcells(schnitznum).ang = []; % JCR 2005-06-23
        for age = 1:schnitzcells(schnitznum).N
            framenum = s.frames(age) - 1; % JCR hack - schnitzedit is 1-based 
            % look up that frame number's index within trackRange
            lincellnumIndex = find(trackRange==framenum);
            cellnum  = s.cellno(age);
            lincellnum{lincellnumIndex}(cellnum) = schnitznum;
        end
    end
    
    if exist(p.flatFieldName)==2,
        load(p.flatFieldName,'zzcim','zzyim',...
            'YFPmolecules_per_YFLunit','CFPmolecules_per_CFLunit');
    else
        disp('No flatfield files!');
    end
    if exist('zzcim')~=1
        disp('zzcim not found, using blank "zzcim".');
        zzcim=[];
    end
    if exist('zzyim')~=1
        disp('zzyim not found, using blank "zzyim".');
        zzyim=[];
    end
    if exist('CFPmolecules_per_CFLunit')==1
        p.CFPmolecules_per_CFLunit = CFPmolecules_per_CFLunit;
    else
        disp('No CFP calibration data - using default=1.');
        p.CFPmolecules_per_CFLunit = 1;
    end
    if exist('YFPmolecules_per_YFLunit')==1
        p.YFPmolecules_per_YFLunit = YFPmolecules_per_YFLunit;
    else
        disp('No YFP calibration data - using default=1.');
        p.YFPmolecules_per_YFLunit = 1;
    end
    %-----------------------added ST 12-12-2005
    if exist('zzgim')~=1
        disp('zzyim not found, using blank "zzgim".');
        zzgim=[];
    end
    if exist('GFPmolecules_per_YFLunit')==1
        p.GFPmolecules_per_GFLunit = GFPmolecules_per_GFLunit;
    else
        disp('No GFP calibration data - using default=1.');
        p.GFPmolecules_per_GFLunit = 1;
    end
    %-----------------------ST
    
    for i = 1:length(trackRange)
        currFrameNum = trackRange(i);
        % i is the frame index into lincellnum and trackRange. 
        % trackRange(i) gives the actual frame number that corresponds 
        % to a given lincellnum{i}.
        % Also, recall that a schnitz' frames array is 1-based per 
        % makeschnitz hack to conform to schnitzedit expectations

        % load segmented image for this frameNum
        name= [p.segmentationDir,p.movieName,'seg',str3(currFrameNum)];
        disp([' Loading frame ', str3(currFrameNum)]);
        clear LNsub Lc;
        creg = [];
        yreg = [];
        greg = [];
        rreg = [];
        exptc = [];
        expty = []; 
        exptg = [];
        exptr = [];
        cbinning = p.camera_bin_size;
        ybinning = p.camera_bin_size;
        gbinning = p.camera_bin_size;
        rbinning = p.camera_bin_size;
        load([name]); % including variables: LNsub, Lc, phsub, timestamp 
        % rect, [and for all colors:] cback, creg, cshift, exptc, gainc
        % keyboard
        if exist('Lc')~=1; Lc=LNsub;end % nitzan added 2005June25
        
        creg=double(creg);yreg=double(yreg);
        if ~isempty(creg)
            creg = normalizeF(creg,cback,zzcim,rect,exptc,p.autoCFL,...
                cbinning,p.CFPmolecules_per_CFLunit);
        end
        if ~isempty(yreg),
          if size(creg,2)==size(yreg,2) & max(creg(:))>0
            yreg = normalizeF(yreg,yback,zzyim,rect,expty,p.yfpOffset,...
                ybinning,p.YFPmolecules_per_YFLunit,creg,p.crosstalkRatioC2Y);
          else % if there's no CFP image we can't subtract crosstalk - assume none.
            yreg = normalizeF(yreg,yback,zzyim,rect,expty,p.yfpOffset,...
                ybinning,p.YFPmolecules_per_YFLunit);
          end
        end
        greg=double(greg);rreg=double(rreg);
        if ~isempty(greg)
            greg = normalizeF(greg,gback,zzgim,rect,exptg,p.autoGFL,...
                gbinning,p.GFPmolecules_per_GFLunit);
        end
        if ~isempty(rreg)
            greg = normalizeF(rreg,rback,zzrim,rect,exptr,p.autoRFL,...
                rbinning,p.RFPmolecules_per_RFLunit);
        end
        rp = regionprops(Lc,'orientation','majoraxislength','minoraxislength','centroid');
               
        % update each schnitz that exists during this frame
        schnitzesForFrame = lincellnum{i};
        nonZeroSchnitzes = schnitzesForFrame;       % nitzan 2005May16
        nonZeroSchnitzes(nonZeroSchnitzes==0)=[];   % nitzan 2005May16
        % nitzan added 2005May24 for the problematic movie #3
        if isfield(p,'absaugOnlyApproved')
            if p.absaugOnlyApproved==1
                nonZeroSchnitzes = nonZeroSchnitzes...
                    (find([(schnitzcells(nonZeroSchnitzes).approved)]));
            end
        end
        % end nitzan added 2005May24 for the problematic movie #3
        for currSchnitzNum = nonZeroSchnitzes       % nitzan 2005May16
        %for currSchnitzNum = schnitzesForFrame     % nitzan 2005May16
            currSchnitz = schnitzcells(currSchnitzNum);
            % figure out index within this schnitz' age-based arrays
            age = find((currSchnitz.frames-1) == currFrameNum);
            if isempty(age)
                error(['lincellnum says schnitz num ' currSchnitzNum ...
                       'exists in frame ' currFrameNum ', but that frame ' ...
                       ' can''t be found in the schnitz' frames array']);
            end
            % get cell number in segmented image for current schnitz & frame
            cellnum = currSchnitz.cellno(age);
            % get cell characteristics from images
            schnitzcells= getcellparams(schnitzcells, rp, Lc, rect, creg,...
                yreg, exptc, expty, currSchnitzNum, age, cellnum); 
            % also add image's datetime timestamp
 %           keyboard
            if (timestamp ~= 'empty')
                schnitzcells(currSchnitzNum).datetime(age) = timestamp;
            else
                movieYear  = str2num(p.movieDate(1:4));
                movieMonth = str2num(p.movieDate(6:7));
                movieDate  = str2num(p.movieDate(9:10));
                fakeMin    = currFrameNum;
                schnitzcells(currSchnitzNum).datetime(age) = ...
                    datenum(movieYear, movieMonth, movieDate, 0, fakeMin, 0);
            end
        end % loop over schnitzesForFrame
    end % loop over frames
else
    disp('quick mode on! skipping fluorescence extraction!');
    % we're in quick mode, so load an existing schnitz file which
    % contains the image-derived fields
    if exist(p.schnitzName)~=2
        error(['Unable to read schnitzcells file ''' schnitzFilePath ...
            ''' which is required for quick mode']);
    end
    load(p.schnitzName);    % 'schnitzcells'
    disp('loading old schnitzcells and updating post-image derived fields...')
end

schnitzcells= add_volume(schnitzcells, p.micronsPerPixel);
% schnitzcells= add_poles(schnitzcells);       % nitzan 2005May16 - this crashed.
% schnitzcells= add_generations(schnitzcells); % nitzan 2005May16 - this crashed.
schnitzcells= add_cyclephase(schnitzcells);
schnitzcells= add_minutes(schnitzcells);

schnitzcells= add_derivs(schnitzcells);
if length([schnitzcells.FCs])
    schnitzcells= add_interps(schnitzcells,'FYsmins','FYs','dFCmins','FY_interps_Cdot');
    schnitzcells= add_interps(schnitzcells,'mins','volume','dFCmins','vol_interps_Cdot');
else
    schnitzcells= add_interps(schnitzcells,'FYsmins','FYs','dFYmins','FY_interps_Ydot');
    schnitzcells= add_interps(schnitzcells,'mins','volume','dFYmins','vol_interps_Ydot');
end
schnitzcells= add_interps(schnitzcells,'mins','phase','dFCmins','phase_interps_Cdot');
schnitzcells= add_interps(schnitzcells,'mins','phase','dFYmins','phase_interps_Ydot');


if isfield(p,'getFig3DemoFields') & p.getFig3DemoFields
    schnitzcells = getFig3DemoFields(schnitzcells,p.micronsPerPixel);
end

save(p.schnitzName, 'schnitzcells');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xreg = SubstrBackGrndF(xreg,xback,zzxim,rect,exptx,Auto_XFL,...
    camera_bin_size,XFPmolecules_per_XFLunit,crosstalkin,crosstalkratio)

xreg=double(xreg)-xback; % in camera_xfp_units; xback also in camera units.
if prod(size(zzxim))>10 
    xreg=xreg./zzxim(rect(1):rect(3), rect(2):rect(4));
end
% exptx=1; %SJT
xreg=xreg/exptx; % in camera_xfp_1sec_units;
xreg=xreg*XFPmolecules_per_XFLunit;
xreg=xreg/(camera_bin_size^2);
xreg=xreg-Auto_XFL; % in units of molecules?
if nargin>8
    xreg = xreg - crosstalkin*crosstalkratio;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xreg = normalizeF(xreg,xback,zzxim,rect,exptx,Auto_XFL,...
    camera_bin_size,XFPmolecules_per_XFLunit,crosstalkin,crosstalkratio)

%xreg=double(xreg)-xback; % in camera_xfp_units; xback also in camera units.
xreg=double(xreg)-0; % SJT in camera_xfp_units; xback also in camera units.

if prod(size(zzxim))>10 
    xreg=xreg./zzxim(rect(1):rect(3), rect(2):rect(4));
end
% exptx=1; %SJT
% xreg=xreg/exptx; % SJT in camera_xfp_1sec_units;
xreg=xreg*XFPmolecules_per_XFLunit;
% xreg=xreg/(camera_bin_size^2); %SJT

xreg=xreg-Auto_XFL; % in units of molecules?

if nargin>8
    xreg = xreg - crosstalkin*crosstalkratio;
end
xreg%AA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ncells= getcellparams(cells, rp, Lc, rect,...
    creg, yreg, exptc, expty, schnitzno, age, cellnum)
ncells= cells;
% find cell
loc= find(Lc == cellnum);
% extract information
ncells(schnitzno).ang(age) = rp(cellnum).Orientation;
ncells(schnitzno).len(age) = rp(cellnum).MajorAxisLength;
ncells(schnitzno).wid(age) = rp(cellnum).MinorAxisLength;
ncells(schnitzno).cenx(age) = rp(cellnum).Centroid(1);
ncells(schnitzno).ceny(age) = rp(cellnum).Centroid(2); 
ncells(schnitzno).cenX(age) = rp(cellnum).Centroid(1) + rect(2) - 1;
ncells(schnitzno).cenY(age) = rp(cellnum).Centroid(2) + rect(1) - 1;    
ncells(schnitzno).SZ(age)= length(loc);
if ~isempty(creg),
    ncells(schnitzno).FC(age)= sum(creg(loc));
    ncells(schnitzno).SC(age)= std(creg(loc));
    ncells(schnitzno).MC(age)= ncells(schnitzno).FC(age)/ncells(schnitzno).SZ(age);
    ncells(schnitzno).medC(age) = median(creg(loc));
else
    ncells(schnitzno).FC(age)= NaN;
    ncells(schnitzno).SC(age)= NaN;
    ncells(schnitzno).MC(age)= NaN;
    ncells(schnitzno).medC(age)= NaN;
end;    
if ~isempty(yreg) | ~isempty(yreg),
    ncells(schnitzno).FY(age)= sum(yreg(loc));
    ncells(schnitzno).SY(age)= std(yreg(loc));
    ncells(schnitzno).MY(age)= ncells(schnitzno).FY(age)/ncells(schnitzno).SZ(age);
    ncells(schnitzno).medY(age) = median(yreg(loc));
else
    ncells(schnitzno).FY(age)= NaN;
    ncells(schnitzno).SY(age)= NaN;
    ncells(schnitzno).MY(age)= NaN;
    ncells(schnitzno).medY(age)= NaN;
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newschnitzcells = add_volume(schnitzcells,microns_per_pixel)

for i = 1:length(schnitzcells),
    lens = schnitzcells(i).len;
    wids = schnitzcells(i).wid;
    rads=wids/2;
    %volume in pixels^3: for a cylinder=pi*(rad^2)*len; 
    %for a cigar shape: , then change to microns^3.
    vols=((4/3)*pi*(rads.^3) + pi*(rads.^2).*(lens-2*rads)); 
    vols=vols*(microns_per_pixel^3);     % in microns^3
    schnitzcells(i).volume  = vols;      % microns cubed
    schnitzcells(i).lengthMicrons = schnitzcells(i).len*microns_per_pixel;
    schnitzcells(i).width_Microns = schnitzcells(i).wid*microns_per_pixel;
end

newschnitzcells = schnitzcells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newschnitz = add_generations(schnitzcells);
for i= 1:length(schnitzcells)
    schnitzcells(i).gen= mygen(schnitzcells, i);
    if min(size(schnitzcells(i).FC))==0,
        schnitzcells(i).gens = [];
    else
        schnitzcells(i).gens= ones(size(schnitzcells(i).FC))*schnitzcells(i).gen;
    end;
end;
newschnitz = schnitzcells;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gen = mygen(schnitzcells, me);
gen = 0;
while (schnitzcells(me).P ~= -1),
    gen = gen+1;
    me = schnitzcells(me).P;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newschnitzcells = add_derivs(schnitzcells);
% modified on 4/22/03 by m.e. & n.r.
% to base the derivaties on time!
% note the complete generality of this and all future functions with
% respect to color.
colors = 'CYGR';
maxcolors = length(colors);
% first do the easy stuff: derivatives within a single schnitz.
for color = 1:maxcolors;
    thisfield = ['F',colors(color)];
    derivname = ['dF',colors(color),'dt']; % fields like dFCdt
    derivolname = ['dF',colors(color),'dvol']; % fields like dFCdvol

    %     derivlenname = ['dF',colors(color),'dlen'];
    %     derivszname = ['dF',colors(color),'dsz'];
    
    derivtimename = ['dF',colors(color),'mins']; % fields like dFCmins (the time of the derivative)
    framename = ['F',colors(color),'sframes'];  % here are some additional fields which contain subsets of the FC, etc. values omitting the NaNs
    timename = ['F',colors(color),'smins'];
    volname = ['F',colors(color),'svols'];
    fluorname = ['F',colors(color),'s'];
    mfluorname = ['M',colors(color),'s'];
    lfluorname = ['L',colors(color),'s'];  
    if isfield(schnitzcells(1),thisfield),
        for i = 1:length(schnitzcells),
            FF = getfield(schnitzcells(i),thisfield);
            f = find(~isnan(FF));
            if isempty(f),
                schnitzcells(i).(derivname)=[];
                schnitzcells(i).(derivolname)=[];
                schnitzcells(i).(framename)=[];
                schnitzcells(i).(fluorname)=[];
                schnitzcells(i).(timename)=[];
                schnitzcells(i).(volname)=[];
                schnitzcells(i).(mfluorname)=[];
                schnitzcells(i).(lfluorname)=[];
                schnitzcells(i).(derivtimename)=[];
            else
                times = schnitzcells(i).mins(f);
                vols = schnitzcells(i).volume(f);
                FFs = FF(f);
                MF = getfield(schnitzcells(i),['M',colors(color)]);
                MFs = MF(f);
                lens = schnitzcells(i).len(f)*12;
                % now we have a vector of the existing times and fluorescence
                % values:
                dtimes = diff(times);
                dvols = diff(vols);
                dFFs = diff(FFs);
                dFFdt = dFFs ./ dtimes;
                dFFdV = dFFs ./ dvols;
                derivtimes = 0.5*(times(1:end-1)+times(2:end));
                % now add these to the schnnitz:
                schnitzcells(i).(derivname)=[dFFdt];
                schnitzcells(i).(derivolname)=[dFFdV];
                schnitzcells(i).(framename)=[schnitzcells(i).frames(f)];
                schnitzcells(i).(fluorname)=[FFs];
                schnitzcells(i).(timename)=[schnitzcells(i).mins(f)];
                schnitzcells(i).(volname)=[schnitzcells(i).volume(f)];
                schnitzcells(i).(mfluorname)=[MFs];
                schnitzcells(i).(lfluorname)=[FFs./lens];
                schnitzcells(i).(derivtimename)=[derivtimes];
            end;
        end;
        for i = 1:length(schnitzcells),
            D = schnitzcells(i).D;
            E = schnitzcells(i).E;
            if (D>0) & (E>0),
                Dfluor = getfield(schnitzcells(D),fluorname);
                Efluor = getfield(schnitzcells(E),fluorname);
                OurFluors = getfield(schnitzcells(i),fluorname);
                Ourtimes = getfield(schnitzcells(i),timename);
                Ourvol = getfield(schnitzcells(i),volname);
                if (length(Dfluor)>0) & (length(Efluor)>0) & (length(OurFluors)>0),
                    TotalFluor = Dfluor(1)+Efluor(1);
                    OurFluor = OurFluors(end);    
                    Ourtime = Ourtimes(end);
                    Ourvol = Ourvol(end);
                    Dtimes = getfield(schnitzcells(D),timename);
                    Dvol = getfield(schnitzcells(D),volname);
                    Evol = getfield(schnitzcells(E),volname);
                    Dtime = Dtimes(1);
                    DEvol = Dvol(1) + Evol(1);
                    delta_time = Dtime - Ourtime;
                    delta_vol = DEvol - Ourvol;
                    deriv_time = (Dtime + Ourtime)/2;
                    deriv = (TotalFluor - OurFluor) / delta_time;
                    derivol = (TotalFluor - OurFluor) / delta_vol;
                    old_deriv = getfield(schnitzcells(i),derivname);
                    old_derivol = getfield(schnitzcells(i),derivolname);
                    old_deriv_time = getfield(schnitzcells(i),derivtimename);
                    old_deriv(end+1) = deriv;
                    old_derivol(end+1) = derivol;
                    old_deriv_time(end+1) = deriv_time;
                    schnitzcells(i).(derivname) = old_deriv;
                    schnitzcells(i).(derivolname) = old_derivol;
                    schnitzcells(i).(derivtimename) = old_deriv_time;
                end;
            end;            
        end;
    end;
end;    
newschnitzcells = schnitzcells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newschnitzcells= add_interps(...
    schnitzcells,oldtimefield,oldfield,interpatfield,newfieldname);
% this adds an interpolated new field, named "newfieldname", which
% interpolates the oldfield values at the times defined by the interpat
% field, which *has to be a time* (in the same units as oldtimefield)
% note that at division points we give up and put a NaN.
%     -----> not any more! as long as there is at least one value
%            for YFP for that cell, we use it throughout under the
%            assumption that it is constant. only then do we divide
%            by the volume which we measure at each time-point anyway.
for i = 1:length(schnitzcells),
    Xi = schnitzcells(i).(interpatfield);
    Xs = schnitzcells(i).(oldtimefield);
    Ys = schnitzcells(i).(oldfield);
    if length(Xs)>0,
        [Xs,Ys] = expand_vector_for_interp1(Xs,Ys,Xi);
%         keyboard;
        schnitzcells(i).(newfieldname) = interp1(Xs,Ys,Xi);;
    else
        schnitzcells(i).(newfieldname) = ...
            NaN * ones(size(Xi));
    end;
    if length(Xi)==0,
        schnitzcells(i).(newfieldname) = [];
    end;
end;
newschnitzcells = schnitzcells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xs2,Ys2] = expand_vector_for_interp1(Xs,Ys,Xi);

if min([0.9999*Xs(1),1.0001*Xs(1),Xi])<Xs(1)
    Xs = [min([0.9999*Xs(1),1.0001*Xs(1),Xi]) , Xs];
    Ys = [Ys(1) , Ys];
end
if max([Xi,0.9999*Xs(end),1.0001*Xs(end)])>Xs(end)
    Xs = [Xs , max([Xi,0.9999*Xs(end),1.0001*Xs(end)])];
    Ys = [Ys , Ys(end)];
end
Xs2=Xs;Ys2=Ys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newschnitzcells = add_cyclephase(schnitzcells);
% adds a cell cycle "phase" variable called phase
% that ranges from 0 on the first frame to 1 on the last frame.
newschnitzcells = schnitzcells;
for i = 1:length(schnitzcells),
    L = length(schnitzcells(i).frames);
    if L==0,
        %disp('zero frames?'); % nitzan 2005May23
    elseif L == 1,
        newschnitzcells(i).phase(1) = 0.5;
    else
        for k = 1:L,
            newschnitzcells(i).phase(k) = (k-1)/(L-1);
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poleschnitz = add_poles(schnitzcells, me, theta0);
% add information about poles and reorder daughters D and E
% D contains the oldest pole
% angle points from old-->new.  i.e. if the old pole was at the origin,
% the new pole would be at (cos(theta), sin(theta))
% first, propagate theta0 until the end of me, 
% then recurse on daughters.
if nargin==1,
    poleschnitz = schnitzcells;
    f = find([schnitzcells.B]==min([schnitzcells.B]));
    for i = 1:length(f),
        poleschnitz = add_poles(poleschnitz,f(i),0);
    end;
    return;
end;
myschnitz = schnitzcells(me);
for i = 1:length(myschnitz.frames),
    mythetas(i) = calctheta(theta0, myschnitz.ang(i));
    theta0 = mythetas(i);
end;
poleschnitz = schnitzcells;
poleschnitz(me).thetas = mythetas;
% now call for our daughters:
if (myschnitz.D > 0) & (myschnitz.E == 0),
    disp('weird -- one daughter only??');
end;
if (myschnitz.D > 0) & (myschnitz.E > 0),
    xD = [poleschnitz(myschnitz.D).cenx(1) poleschnitz(myschnitz.D).ceny(1)];
    xE = [poleschnitz(myschnitz.E).cenx(1) poleschnitz(myschnitz.E).ceny(1)];
    % first consider cell D:
    % vec from D to E:
    vec = (xE-xD)/sqrt(sum((xE-xD).^2));
    % two possible angles:
    [theta1, theta2] = ang2theta(poleschnitz(myschnitz.D).ang(1));
    vec1 = [cos(theta1) -sin(theta1)];
    vec2 = [cos(theta2) -sin(theta2)];
    dot1 = acos(sum(vec1.*vec));
    dot2 = acos(sum(vec2.*vec));
    if dot1 < dot2
        thetaD= theta1;
    else
        thetaD= theta2;
    end;
    % second consider cell E:
    % vec from E to D:
    vec = (xD-xE)/sqrt(sum((xE-xD).^2));
    % two possible angles:
    [theta1, theta2] = ang2theta(poleschnitz(myschnitz.E).ang(1));
    vec1 = [cos(theta1) -sin(theta1)];
    vec2 = [cos(theta2) -sin(theta2)];
    dot1 = acos(sum(vec1.*vec));
    dot2 = acos(sum(vec2.*vec));
    if dot1 < dot2
        thetaE= theta1;
    else
        thetaE= theta2;
    end;
    % order daughters dependent on mother cell
    %     mdang= abs([theta0 - thetaE, theta0 - thetaD]);
    deltaD = deltatheta(theta0,thetaD);
    deltaE = deltatheta(theta0,thetaE);
    daughterD= myschnitz.D;
    daughterE= myschnitz.E;
    savethetaD = thetaD;
    savethetaE = thetaE;
    if deltaE < deltaD,
        % myschnitz.D contains the oldest pole
        poleschnitz(me).D = daughterE;
        poleschnitz(me).E = daughterD;
        thetaD = savethetaE;
        thetaE = savethetaD;
    end;
    poleschnitz = add_poles(poleschnitz, poleschnitz(me).D, thetaD);  
    poleschnitz = add_poles(poleschnitz, poleschnitz(me).E, thetaE); 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newtheta = calctheta(oldtheta, curang); % used by add_poles
% curang is from regionprops "orientation", in degrees, which goes 
% from -90-->+90
% theta is from 0-->2pi
[curangrad1, curangrad2] = ang2theta(curang);
delta1 = deltatheta(oldtheta,curangrad1);
delta2 = deltatheta(oldtheta,curangrad2);
if (delta1) < (delta2),
    newtheta = curangrad1;
else
    newtheta = curangrad2;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta1, theta2] = ang2theta(ang); % used by add_poles
curangrad = pi * ang/180;
% curangrad is from 0-->2pi
if curangrad < 0
    curangrad1= curangrad + 2*pi;
    curangrad2= curangrad - pi + 2*pi;
else
    curangrad1 = curangrad;
    curangrad2 = curangrad + pi;
end;
theta1 = curangrad1;
theta2 = curangrad2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = deltatheta(theta1,theta2); % used by add_poles
delta = min(abs(theta1 - [theta2 theta2-2*pi theta2+2*pi]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function schnitzcells2 = add_minutes(schnitzcells);
s = sort([schnitzcells.datetime]);
firstmoment = s(1);
for i = 1:length(schnitzcells),
    tdiff = schnitzcells(i).datetime - firstmoment;
    [y,m,d,h,mi,s] = datevec(tdiff);    
    schnitzcells(i).mins = s/60 + mi + h*60 + d*24*60 + m *30*24*60+y*365*24*60;
end;
schnitzcells2 = schnitzcells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function schnitzcells = getFig3DemoFields(schnitzcells,microns_per_pixel);
maxYvalM = max([schnitzcells.MYs]);
maxCvalM = max([schnitzcells.MCs]);
maxYval = max([schnitzcells.FYs]);
maxCval = max([schnitzcells.FCs]);
for i=1:length(schnitzcells);
    schnitzcells(i).demo_FYsNorm = schnitzcells(i).FYs/maxYval;
    schnitzcells(i).demo_FCsNorm = schnitzcells(i).FCs/maxCval;
    schnitzcells(i).demo_MYsNorm = schnitzcells(i).MYs/maxYvalM;
    schnitzcells(i).demo_MCsNorm = schnitzcells(i).FCs/maxCvalM;
    if length(schnitzcells(i).FCs)>3
        schnitzcells(i).demo_dFCvals = [diff(schnitzcells(i).demo_FCsNorm),NaN];
        schnitzcells(i).demo_dFCtime = [(schnitzcells(i).FCsmins(1:end-1)+...
            schnitzcells(i).FCsmins(2:end))/2,NaN];
    else
        schnitzcells(i).demo_dFCvals = [NaN];
        schnitzcells(i).demo_dFCtime = [NaN];
    end
    schnitzcells(i).demo_dFCvals(schnitzcells(i).demo_dFCvals<=0)=NaN;
    if length(schnitzcells(i).FYs)>1
        schnitzcells(i).demo_dFYvals = 5*[diff(schnitzcells(i).demo_FYsNorm),NaN];
        schnitzcells(i).demo_dFYtime = [(schnitzcells(i).FYsmins(1:end-1)+...
            schnitzcells(i).FYsmins(2:end))/2,NaN];
    else
        schnitzcells(i).demo_dFYvals = [NaN];
        schnitzcells(i).demo_dFYtime = [NaN];
    end
    schnitzcells(i).demo_dFYvals(schnitzcells(i).demo_dFYvals<=0)=NaN;
end