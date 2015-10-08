
%% INSTRUCTIONS

% Set the filepath of the configuration file (by setting CONFIGFILEPATH 
% below), and then run this script by pressing "run and advance" for each 
% section.

%% It is very important to supply the path to the correct dataset here!
CONFIGFILEPATH = 'F:\A_Tans2_step1_incoming_not_backed_up\2015-09-16\Schnitzcells_config_fullAnalysis_2015_09_16.xlsx';

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
    if strcmp(p.(currentFluor),'None')
        disp([currentFluor ' not set, skipping']);
        break
    end
    
    % Commands that need to be executed per color
    % ===
    
    % Load background (mistakenly called flatfield; legacy), shading and
    % replace matrices
    load(settings.fluorCorrectionImagePaths{colorIdx}, 'flatfield', 'shading', 'replace');
    
    % Finding shifts.
    optimalShift = DJK_getFluorShift_anycolor(p,'manualRange', settings.frameRangeFull,'fluorcolor',currentFluor,'maxShift',MAXSHIFT);
    % Correct images (shading, background).
    DJK_correctFluorImage_anycolor(p, flatfield, shading, replace,'manualRange', settings.frameRangeFull,  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF),'fluorcolor',currentFluor,'minimalMode',0);

end
    
% Let user hear we're done
mysound=load('gong'); sound(mysound.y);


%%











