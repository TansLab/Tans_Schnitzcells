
%% INSTRUCTIONS

% Set the filepath of the configuration file (by setting CONFIGFILEPATH 
% below), and then run this script by pressing "run and advance" for each 
% section.

%% It is very important to supply the path to the correct dataset here!
CONFIGFILEPATH = 'F:\A_Tans2_step1_incoming_not_backed_up\2015-09-16\Schnitzcells_config_fullAnalysis_2015_09_16_pos1.xlsx';

%% Read configuration settings from excel file.

% One can execute this section again to reload settings 
% Note that these are then not immediate parsed to the "p" struct.
settings = MW_readsettingsfromexcelfile(CONFIGFILEPATH)

%% Double check whether to continue
% This step is mainly to prevent accidentally running whole script 
myAnswer = questdlg(['Loaded ' CONFIGFILEPATH ', do you want to continue?'],'Confirmation required.','Yes','No','No');
if strcmp(myAnswer, 'No') || strcmp(myAnswer, 'Cancel')
    error('Analysis aborted.');
end

%% Now make a vector with parameter values that schnitzcells scripts can handle.

disp('Now creating ''p'' struct from settings struct.');

% Create the p vector which holds all parameters and is fed into, and also
% returned by most functions of the schnitzcells analysis software.
% Use "settings" struct as a base.
% DJK_initschnitz only checks for errors and adds a few parameters based on
% the already given paramteres.
% ===
% TODO this can be done more elegantly (but note that existence of two
% vectors, "settings" and "p", allows user to update "settings" vector 
% intermediately.
p = DJK_initschnitz(settings.positionName,settings.movieDate,'e.coli.amolf','rootDir',...
    settings.rootDir, 'cropLeftTop',settings.cropLeftTop, 'cropRightBottom',settings.cropRightBottom,...
    'fluor1',settings.fluor1,'fluor2',settings.fluor2,'fluor3',settings.fluor3,...
    'setup',settings.setup,'softwarePackage',settings.softwarePackage,'camera',settings.camera)

% Manually make sure image dir is correct
% (This is done to accomodate cropping.)
p.imageDir = [settings.rootDir settings.movieDate '\' settings.positionName '\']

%% Crop images

% Crop images 
% =========================================================================

% Crop images
% Puts cropped images in new directory, which gets the suffix 
DJK_cropImages_3colors(p, settings.frameRangeFull, settings.cropLeftTop, ...
    settings.cropRightBottom, 'cropName', [settings.positionName settings.cropSuffix]);

% Update p accordingly.
% =========================================================================
% Load this setting to continue an analysis from here
% =========================================================================
% (After loading settings from excel file.)
p = DJK_initschnitz([settings.positionName settings.cropSuffix],settings.movieDate,'e.coli.amolf','rootDir',...
    settings.rootDir, 'cropLeftTop',settings.cropLeftTop, 'cropRightBottom',settings.cropRightBottom,...
    'fluor1',settings.fluor1,'fluor2',settings.fluor2,'fluor3',settings.fluor3,...
    'setup',settings.setup,'softwarePackage',settings.softwarePackage,'camera',settings.camera)
% =========================================================================

disp('Done cropping');

%% Segmenation
% =========================================================================
% ADVANCED PARAMETER SETTINGS that can be useful:
% p.useFullImage=1; % forces to use full image for segmentation
% p.overwrite=0; % overwrite existing files

PN_segmoviephase_3colors(p,'segRange', settings.frameRangeFull,'slices', settings.slices,...
    'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
    settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
    'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

PN_copySegFiles(p,'segRange', settings.frameRangeFull,'slices', settings.slices,...
    'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
    settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
    'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

% ADVANCED PARAMETER
p.showAll = 1; % show all segmented frames to user

% Let user now section is done by sound
mysound=load('gong'); sound(mysound.y);

% Start the manual checking by user
PN_manualcheckseg(p,'manualRange',settings.frameRangeFull,'override',0,'assistedCorrection',1); % ADVANCED PARAMETERS in fn arguments


%% Perform quick analysis of segmentation
DJK_analyzeSeg(p,'manualRange',settings.frameRangeFull,'onscreen',1);

%% Perform tracking, check tracking, make problem movie and correct
% Iterate this section until no problems are listed any more in:
% \analysis\tracking\manualRangeX_Y\posZcrop-tracking.txt
% iterate these commands (also partially redundant w. above) to correct segmentation
p.overwrite=0; % ADVANCED SETTING
p.showAll = 1; % show all segmented frames to user

% Alternative trackers one can try:
% slow but more robust
%NW_tracker_centroid_vs_area(p,'manualRange', [1:244]); 
% fast but fails often
%MW_tracker(p,'manualRange', [1:244]); 
% Default one
DJK_tracker_djk(p,'manualRange', settings.frameRangeFull); % default tracker

% Find problem cells
problems = DJK_analyzeTracking(p,'manualRange', settings.frameRangeFull, 'pixelsMoveDef', 15, 'pixelsLenDef', [-4 13]);
p.problemCells = problems;

% Create lookup table for frame label <-> schnitz label
p.slookup=MW_makeslookup(p);

% Let user know we're done
mysound=load('gong'); sound(mysound.y);

% Manual check again 
% Since "p" now contains the lookup table, and problemcells is defined, it
% will highlight the problemcells.
% Note that when one performs manual corrections to a frame, the mapping is
% not correct any more.
PN_manualcheckseg(p,'manualRange',settings.frameRangeFull,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted

%% Once above section has been iterated to satisfaction make sure to remove problemCells field from p
% (Because movies later also work with this.)
p =rmfield(p,'problemCells');

%% Make movie that allows you to double check the segmentation & tracking

% User info
if isfield(p, 'problemCells')
    disp('p.problemCells field exists, will make movie of those frames only.');
end
    
% Make movie
DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1,'manualRange',settings.frameRangeFull);

% User info
disp(['Result written to ' p.analysisDir 'movies\']);

%{
% Other movie making options
DJK_makeMovie(p, 'tree', 'schAll', 'stabilize', 1,'problemCells',problems);
DJK_makeMovie (p, 'tree', 'cellno', 'stabilize', 1,'problemCells',problems);
p=rmfield(p,'problemCells')
DJK_makeMovie (p, 'tree', 'cellno', 'stabilize', 1);
%DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1);
%}

mysound=load('gong'); sound(mysound.y);

%% Create backup of segtrack
NW_Backup_SegTrack(p);

%% Start of correcting fluor
MAXSHIFT = 10;

% Initialization, applies to all colors
% ===
NW_initializeFluorData(p,'manualRange', settings.frameRangeFull);
% Load PSF (color-independent)
load(settings.fluorPointSpreadFunctionPath, 'PSF');
    
% Loop over all colors
% ===
fluorColors = {'fluor1','fluor2','fluor3'}; % ugly but compatible w. legacy - MW
for colorIdx = 1:3
    
    % Select current string w. color identifier (fluor1, ..)
    currentFluor = fluorColors{colorIdx};
    if strcmp(p.(currentFluor),'none') || strcmp(p.(currentFluor),'None')
        disp([currentFluor ' not set, assuming you don''t have more fluor colors.']);
        break
    end
    
    % Commands that need to be executed per color
    % ===
    
    % Load background (mistakenly called flatfield; legacy), shading and
    % replace matrices
    load(settings.fluorCorrectionImagePaths{colorIdx}, 'flatfield', 'shading', 'replace');
    
    % Finding shifts.
    disp('Looking for shifts');
    optimalShift = DJK_getFluorShift_anycolor(p,'manualRange', settings.frameRangeFull,'fluorcolor',currentFluor,'maxShift',MAXSHIFT);
    disp('Correcting');
    % Correct images (shading, background).
    DJK_correctFluorImage_anycolor(p, flatfield, shading, replace,'manualRange', settings.frameRangeFull,  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF),'fluorcolor',currentFluor,'minimalMode',0);
    
end

% Add fluor information to schnitz struct (schnitzcells struct is saved).
disp('Adding fluor info to schnitzcells structure.')
DJK_compileSchnitzImproved_3colors(p,'quickMode',0);

% Let user hear we're done
mysound=load('gong'); sound(mysound.y);

%% Now calculate lengths and add them to schnitzcells struct
% Some of the fluorescence analysis also needs the fields generated here

% ADVANCED PARAMETERS
%p.schnitzNum='all';
%p.onscreen=1;

DJK_addToSchnitzes_length(p);

%8 Add correct mu
DJK_addToSchnitzes_mu(p, 'frameSizes', settings.muWindow);
%DJK_addToSchnitzes_mu(p, 'onScreen', 1, 'frameSizes', [5, 9]);

mysound=load('gong'); sound(mysound.y);

%% Continue with fluor again now

% Again loop over colors
fluorColors = {'fluor1','fluor2','fluor3'}; % ugly but compatible w. legacy - MW
for colorIdx = 1:3
    
    % Adminstration
    % Select current string w. color identifier (fluor1, ..)    
    currentFluor = fluorColors{colorIdx};
    if strcmp(p.(currentFluor),'none') || strcmp(p.(currentFluor),'None')
        disp([currentFluor ' not set, assuming you don''t have more fluor colors.']);
        break
    end
    
    % Actual code that needs to be executed
    % ===
    % Add fluor values to schnitcells struct
    DJK_addToSchnitzes_fluor_anycolor(p, 'onScreen', 0,'colorNormalize',[0 300], 'fluorcolor',currentFluor,'minimalMode',0);

    % Add old style rate vector (not used any more)
    schnitzcells = DJK_addToSchnitzes_fluorRate_phase(p,p.(currentFluor),'5');
    
    % Add new style rate vectors.
    schnitzcells = MW_fluorRate_anycolor(schnitzcells,p.(currentFluor),'length_fitNew');    
    
    % Note PN_smooth_field also adds new field (so you always get
    % smoothing? TODO?! MW)
    myFieldName = ['d' upper(p.(currentFluor)) '5_sum_dt'];
    schnitzcells = PN_smooth_field(schnitzcells,myFieldName,'extensive','length_fitNew');
    NW_saveSchnitzcells(p,schnitzcells); % could be integrated to PN_smooth_field
        
    % only necessary 1 time
    if strcmp(currentFluor, fluorColors{1}) 
        schnitzcells = DJK_addToSchnitzes_predictedValues(schnitzcells, 'phase', 'length_fitNew', 'phase2', [0 1]);
        NW_saveSchnitzcells(p,schnitzcells); % could be integrated to DJK_addToSchnitzes_predictedValues
    end
       
    % Add some more fields
    % Char coding for current color
    X = upper(p.(currentFluor));  
    schnitzcells = MW_addToSchnitzes_atXandDX(schnitzcells, 'phase2',X);
    schnitzcells = MW_addToSchnitzes_atXandDX(schnitzcells, 'time',X);
    schnitzcells = PN_addToSchnitzes_Phase_at_TimeField(schnitzcells,...
        ['phase2'],...
        ['d' X '5_time']); % dX5_time
    schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
        ['d' X '5'],... % dX5
        ['phase2_at_d' X '5_time']);
    schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
        [X '6_mean'],... % X6mean
        ['phase2_at' X '']); % etc..
    
    for muWindowSize = settings.muWindow
        currentMuWindowSizeStr = num2str(muWindowSize);
        schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
            ['muP' currentMuWindowSizeStr '_fitNew'],...
            ['phase2_at' X '']);
    end
    
    schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
        ['d' X '5_sum_dt_s'],...
        ['phase2_atd' X '']);
    schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
        ['d' X '5_sum_dt'],...
        ['phase2_atd' X '']);
    schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
        [X '5_mean'],...
        ['phase2_at' X '']);

    % Cropping vectors
    for muWindowSize = settings.muWindow
        
        currentMuWindowSizeStr = num2str(muWindowSize);
        
        schnitzcells=PN_addToSchnitzes_Phase_at_TimeField_Mu(schnitzcells,...
            ['muP' currentMuWindowSizeStr '_fitNew_all'],...
            ['d' X '5_time'],...
            ['muP' currentMuWindowSizeStr '_fitNew_atd' X '5']);
        schnitzcells = DJK_addToSchnitzes_cycleCor( schnitzcells,...
            ['muP' currentMuWindowSizeStr '_fitNew_atd' X '5'],...
            ['phase2_at_d' X '5_time']);
        schnitzcells=PN_addToSchnitzes_Phase_at_TimeField_Mu(schnitzcells,...
            ['muP' currentMuWindowSizeStr '_fitNew_all'],...
            ['time_atd' X ''],...
            ['muP' currentMuWindowSizeStr '_fitNew_atd' X '5sumdt']);
        schnitzcells = DJK_addToSchnitzes_cycleCor( schnitzcells,...
            ['muP' currentMuWindowSizeStr '_fitNew_atd' X '5sumdt'],...
            ['phase2_atd' X '']);

    end
    
    % Save updated schnitzcells structure
    NW_saveSchnitzcells(p,schnitzcells); 
end

disp('Done with fluor part II.');

%% DONE

% We are done with settings up schnitzcells structure! All data is saved.
% You can proceed with analyzing the data..

% Load this dataset using
% =========================================================================
p = DJK_initschnitz([settings.positionName settings.cropSuffix], settings.movieDate,'e.coli.AMOLF','rootDir',...
 settings.rootDir, 'cropLeftTop', settings.cropLeftTop, 'cropRightBottom', settings.cropRightBottom,...
     'fluor1',settings.fluor1,...
     'fluor2',settings.fluor2,...
     'fluor3',settings.fluor3,...
     'setup',settings.setup,...
     'softwarePackage',settings.softwarePackage,...
     'camera',settings.camera);
     
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
% =========================================================================



