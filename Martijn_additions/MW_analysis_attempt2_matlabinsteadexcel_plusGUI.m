

%% INSTRUCTIONS

% Set the filepath of the configuration file (by setting CONFIGFILEPATH 
% below), and then run this script by pressing "run and advance" for each 
% section.

%% TODOs

% (1)
% Note sure whether preliminary analysis should save as much to the
% schnitzcells .mat file as it does now. Maybe it is smarter to separate
% the info that analyis saves to a separate file. (This might have been
% implemented better in the excell files with command lists that were
% previously used..)

%% It is very important to supply the path to the correct dataset here!
DATAPATH = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-12-02\';
CONFIGFILENAME = 'Schnitzcells_Analysis_Config_2015_12_02_pos2col1.xlsx';
ANALYSISTYPE = 2; % 1 = preliminary, 2 = full
LOADDATASETATEND = 1;

%% Parameters you SHOULD NOT change
EXCELREADSTART = 14; % line where list of parameters starts in Excel file.

%% Read configuration settings from excel file.

% One can execute this section again to reload settings 
% Note that these are then not immediate parsed to the "p" struct.
[settings, alldata] = MW_readsettingsfromexcelfile([DATAPATH CONFIGFILENAME])


%% Double check whether to continue
% This step is mainly to prevent accidentally running whole script 
analysisTypes = {'PRELIMINARY','FULL'}

myAnswer = questdlg(['Loaded ' [DATAPATH CONFIGFILENAME] ', for ' analysisTypes{ANALYSISTYPE} ' analysis do you want to continue?'],'Confirmation required.','Yes','No','No');
if strcmp(myAnswer, 'No') || strcmp(myAnswer, 'Cancel')
    error('Analysis aborted.');
end

%% Now make a vector with parameter values that schnitzcells scripts can handle. (and some misc. other admin)

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

% Set framerange according to analysis type
if ANALYSISTYPE == 1 % fast analysis
    currentFrameRange = settings.frameRangePreliminary;
elseif ANALYSISTYPE == 2 % full analysis
    currentFrameRange = settings.frameRangeFull;
else
    error('No analysis type speficied');
end

%% Some more hard coded parameters

% Output directory for figures etc.
MYOUTPUTDIR = [p.dateDir 'outputSummary\'];

%% Crop images

% Crop images 
% =========================================================================

% Crop images
% Puts cropped images in new directory, which gets the suffix 
if ANALYSISTYPE == 1
    % Determine crop area
     [selectedLeftTop,selectedRightBottom] = MW_determinecroparea(p, settings.frameRangePreliminary);

    % Close figure from crop popup
    close(gcf);
     
    % set new croparea and save crop area to excel file
    myAnswer = questdlg(['Start cropping? And save selection to Excel file (close it first)?'],'Confirmation required.','Save,use,crop','Crop using old','Abort!','Save,use,crop');
    if strcmp(myAnswer, 'Save,use,crop')
        settings.cropLeftTop = selectedLeftTop;
        settings.cropRightBottom = selectedRightBottom;
        
        ExcelcropLeftTopIndex = find(strcmp({alldata{:,1}},'cropLeftTop'))+EXCELREADSTART-1; % find line w. cropLeftTop field.
        xlswrite([DATAPATH CONFIGFILENAME],{mat2str(settings.cropLeftTop)},['B' num2str(ExcelcropLeftTopIndex) ':B' num2str(ExcelcropLeftTopIndex) '']); % write value to it
        ExcelcropRightBottomIndex = find(strcmp({alldata{:,1}},'cropRightBottom'))+EXCELREADSTART-1; % find line w. cropRightBottom field.    
        xlswrite([DATAPATH CONFIGFILENAME],{mat2str(settings.cropRightBottom)},['B' num2str(ExcelcropRightBottomIndex) ':B' num2str(ExcelcropRightBottomIndex) '']); % write value to it
    end
    
    % cropping itself
    if strcmp(myAnswer, 'Save,use,crop') | strcmp(myAnswer, 'Crop using old')
        % Crop images
        DJK_cropImages_3colors(p, settings.frameRangePreliminary, settings.cropLeftTop, ...
            settings.cropRightBottom, 'cropName', [settings.positionName settings.cropSuffix]);    
    end
    
elseif ANALYSISTYPE == 2    
    DJK_cropImages_3colors(p, settings.frameRangeFull, settings.cropLeftTop, ...
        settings.cropRightBottom, 'cropName', [settings.positionName settings.cropSuffix]);
else
    error('Analysis type not specified');
end

disp('Done cropping');

%% Update p accordingly.
% =========================================================================
% Load this setting later if you want to skip cropping
% =========================================================================
% (After loading settings from excel file.)
p = DJK_initschnitz([settings.positionName settings.cropSuffix],settings.movieDate,'e.coli.amolf','rootDir',...
    settings.rootDir, 'cropLeftTop',settings.cropLeftTop, 'cropRightBottom',settings.cropRightBottom,...
    'fluor1',settings.fluor1,'fluor2',settings.fluor2,'fluor3',settings.fluor3,...
    'setup',settings.setup,'softwarePackage',settings.softwarePackage,'camera',settings.camera)
% =========================================================================

disp('p file loaded');

%% Segmenation
% =========================================================================
% ADVANCED PARAMETER SETTINGS that can be useful:
% p.useFullImage=1; % forces to use full image for segmentation
% p.overwrite=1; % overwrite existing files

PN_segmoviephase_3colors(p,'segRange', currentFrameRange,'slices', settings.slices,...
    'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
    settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
    'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

PN_copySegFiles(p,'segRange', currentFrameRange,'slices', settings.slices,...
    'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
    settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
    'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

% ADVANCED PARAMETER
p.showAll = 1; % show all segmented frames to user

% Let user now section is done by sound
mysound=load('gong'); sound(mysound.y);

% To manually redo 1 frame, execute this code manually:
if 0    
    p.overwrite=1;
    TOREDOFRAME = 52:156;
    theOriginalFrameRange = currentFrameRange; currentFrameRange = TOREDOFRAME;
    
    SLICESTEMPORARY = [1 2 3]; % default [1 2 3] instead of settings.slices
    LOGSMOOTHINGTEMPORARY = 2; % default 2; instead of settings.LoG_Smoothing
    MINDEPTHTEMP = 5;% default 5;
    MINCELLAREA = 250; % default 250
    
    PN_segmoviephase_3colors(p,'segRange', currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', MINDEPTHTEMP,'neckDepth', settings.neckDepth);

    PN_copySegFiles(p,'segRange', currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', MINDEPTHTEMP,'neckDepth', settings.neckDepth);    
   
    currentFrameRange = theOriginalFrameRange;
    p.overwrite=0;
end

disp('Done with segmentation.');

%% Start the manual checking by user
if ANALYSISTYPE==1
    assistedYesNo=0;
elseif ANALYSISTYPE==2
    assistedYesNo=1;
end

PN_manualcheckseg(p,'manualRange',currentFrameRange,'override',0,'assistedCorrection',assistedYesNo); % ADVANCED PARAMETERS in fn arguments


%% Perform quick analysis of segmentation
DJK_analyzeSeg(p,'manualRange',currentFrameRange,'onscreen',1,'DJK_saveDir',MYOUTPUTDIR);
%close(gcf);

%% Perform tracking, check tracking, make problem movie and correct
% Iterate this section until no problems are listed any more in:
% \analysis\tracking\manualRangeX_Y\posZcrop-tracking.txt
% iterate these commands (also partially redundant w. above) to correct segmentation

if ANALYSISTYPE == 1 % fast
    DJK_trackcomplete(p,'trackRange',currentFrameRange,'trackMethod','singleCell');
elseif ANALYSISTYPE == 2 % full
    p.overwrite=0; % ADVANCED SETTING
    p.showAll = 1; % show all segmented frames to user

    % Alternative trackers one can try
    % This is quite useful if errors occur. manualrange can be edited to
    % two frames only. Set p.overwrite to 1. Then re-run this whole 
    % section to update general tracking file. (p.overwrite is reset
    % automatically above).
    % ===
    % slow but more robust:
    %NW_tracker_centroid_vs_area(p,'manualRange', [1:244]); 
    % fast but fails often:
    %MW_tracker(p,'manualRange', [1:244]); 
    % If all else fails:
    % edit MW_helperforlinkingframes
    % Default one:
    DJK_tracker_djk(p,'manualRange', currentFrameRange); % default tracker

    % Find problem cells
    problems = DJK_analyzeTracking(p,'manualRange', currentFrameRange, 'pixelsMoveDef', 15, 'pixelsLenDef', [-4 13]);
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
    PN_manualcheckseg(p,'manualRange',currentFrameRange,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted
else
    error('ANALYSISTYPE not speficied');
end

%% Once above section has been iterated to satisfaction make sure to remove problemCells field from p
% (Because movies later also work with this.)
if ANALYSISTYPE == 2
    p =rmfield(p,'problemCells');
end

%% Make movie that allows you to double check the segmentation & tracking

if ANALYSISTYPE == 2

    % User info
    if isfield(p, 'problemCells')
        disp('p.problemCells field exists, will make movie of those frames only.');
    end

    % Make movie
    DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1,'manualRange',currentFrameRange);

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
end

%% Create backup of segtrack
if ANALYSISTYPE == 2
    NW_Backup_SegTrack(p);
end

%% Fluor part I: Start of correcting fluor
MAXSHIFT = 10;

% Initialization, applies to all colors
% ===
NW_initializeFluorData(p,'manualRange', currentFrameRange);
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
    optimalShift = DJK_getFluorShift_anycolor(p,'manualRange', currentFrameRange,'fluorcolor',currentFluor,'maxShift',MAXSHIFT);
    disp('Correcting');
    % Correct images (shading, background).
    DJK_correctFluorImage_anycolor(p, flatfield, shading, replace,'manualRange', currentFrameRange,  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF),'fluorcolor',currentFluor,'minimalMode',0);
    
end

% Add fluor information to schnitz struct (schnitzcells struct is saved).
disp('Adding fluor info to schnitzcells structure.')
DJK_compileSchnitzImproved_3colors(p,'quickMode',0);

% Let user hear we're done
mysound=load('gong'); sound(mysound.y);

disp('Done correcting fluor part I.');

%% Lengths: Now calculate lengths and add them to schnitzcells struct
% Some of the fluorescence analysis also needs the fields generated here

% ADVANCED PARAMETERS
%p.schnitzNum='all';
%p.onscreen=1;

DJK_addToSchnitzes_length(p);

%8 Add correct mu
DJK_addToSchnitzes_mu(p, 'frameSizes', settings.muWindow);
%DJK_addToSchnitzes_mu(p, 'onScreen', 1, 'frameSizes', [5, 9]);

% Let user know we're done
mysound=load('gong'); sound(mysound.y);
disp('Done calculating lengths.');

%% Fluor part II: Continue with fluor again now

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

% Let user know we're done
mysound=load('gong'); sound(mysound.y);
disp('Done with fluor part II.');

%% For preliminary analysis, create some data output
%  ===

if ANALYSISTYPE==1 % fast

    % Obtain schnitzcells
    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

    % Fit mu
    [fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'onScreen', 1,'fitTime',[0 10000],'DJK_saveDir',MYOUTPUTDIR);
    close(gcf);
    %[fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 900], 'onScreen', 1,'fitTime',[0 10000]);

    % Again loop over colors, plot signal for each
    fluorColors = {'fluor1','fluor2','fluor3'}; % ugly but compatible w. legacy - MW
    fitFluorMean = nan(1,3); fitFluorVariance = nan(1,3);
    for colorIdx = 1:3

        % Adminstration
        % Select current string w. color identifier (fluor1, ..)    
        currentFluor = fluorColors{colorIdx};
        if strcmp(p.(currentFluor),'none') || strcmp(p.(currentFluor),'None')
            disp([currentFluor ' not set, assuming you don''t have more fluor colors.']);
            break
        end

        % plot fluor behavior
        fluorFieldName = [upper(p.(currentFluor)(1)) '5_mean_all'];
        [fitFluorMean(colorIdx), fitFluorVariance(colorIdx)] = DJK_plot_avColonyOverTime(p, schnitzcells, fluorFieldName, 'fitTime', fitTime, 'onScreen', 1,'DJK_saveDir',MYOUTPUTDIR);
        close(gcf);
        %DJK_plot_avColonyOverTime(p, schnitzcells, 'C5_mean_all', 'xlim', [0 900], 'ylim', [0 150], 'fitTime', fitTime, 'onScreen',1);
    end

    numberofSchnitzesInLast = MW_countcellsinlastframe(schnitzcells);

    % Write to summary file
    summaryParametersNames = {'#Schnitzes in last frame', 'fitMu', 'fitTime(1)', 'fitTime(2)', 'fitFluorMean1', 'fitFluorVariance1','fitFluorMean2', 'fitFluorVariance2','fitFluorMean3', 'fitFluorVariance3'}
    summaryParameters = [numberofSchnitzesInLast, fitMu, fitTime, fitFluorMean, fitFluorVariance]

    % obtain the number of this position (e.g. pos1 => 1)
    posNumber = str2num(p.movieName(regexp(p.movieName,'[*\d]')))
    
    lineToWriteTo = num2str(posNumber+1);
    
    % Output to excel sheet
    xlswrite([MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{'Identifier'},['A1:A1'])
    xlswrite([MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{[p.movieDate '_' p.movieName]},['A' lineToWriteTo ':' 'A' lineToWriteTo])
    
    for i = 1:numel(summaryParameters)
        letterToWriteTo = ['' i+64+1]; % +1 since start at B
        %disp(['writing ' [letterToWriteTo lineToWriteTo]]);
        xlswrite([MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{summaryParametersNames{i}},[letterToWriteTo '1:' letterToWriteTo '1'])
        xlswrite([MYOUTPUTDIR 'summaryParametersPreliminary.xls'],summaryParameters(i),[letterToWriteTo lineToWriteTo ':' letterToWriteTo lineToWriteTo])
    end

    % Save output to .mat file (update the matfile)
    if exist([MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'file') == 2
        load([MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata');
    end
    thedata(posNumber).summaryParameters = summaryParameters;
    thedata(posNumber).settings = settings;
    thedata(posNumber).p = p; % "settings" and "p" are a bit redundant for historic reasons.
    save([MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata','summaryParametersNames');
    
disp('Done making summary preliminary analysis.');
disp('Check out MW_summaryplotspreliminaryanalysis.m for more plots when all positions done');
    
end


%% DONE / LOADING DATA

% We are done with settings up schnitzcells structure! All data is saved.
% You can proceed with analyzing the data..

% To load after full analysis is done:
if exist('LOADDATASETATEND','var'), if LOADDATASETATEND

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
    
end, end

%% More analyses for full analysis
PLOTSCATTER=0; % activate to plot all scatter plots

if ANALYSISTYPE==2 % full
    
    % Settings up some parameters
    p.myID = settings.myID    
    myOutputFolder = [DATAPATH  '\' p. movieDate  '_' p.movieName '_' p.myID  '\'];
    
    p.NW_saveDir = [myOutputFolder 'misc\'];  % To send additional output to
    p.DJK_saveDir = [myOutputFolder 'misc\']; % To send additional output to
    
    % Location of .mat file containing schnitzcells struct
    myDataFile = [DATAPATH '\' p.movieName  '\data\' p.movieName '-Schnitz.mat'];
   
    % Load datafile
    load(myDataFile);

    % Some more parameter renaming
    myFile = myDataFile; EXPORTFOLDER = myOutputFolder; FITTIME = settings.fitTimeMu;    
    
    % Set up export directory
    if ~exist(EXPORTFOLDER,'dir'), mkdir(EXPORTFOLDER); end
    
    % Some plots of cell growth
    % ===
    % execute script
    MW_growthplots
    
    % create correlation functions R(concentration, growth)
    % ===
    
    % Set up appropriate field names for R(concentration,mu)
    associatedFieldNames = {settings.timeFieldName, settings.fluorFieldName, settings.muFieldName};
    % obtain some settings from Excel file
    badSchnitzes = settings.badSchnitzes; alreadyRemovedInMatFile = settings.alreadyRemovedInMatFile;
    myID = settings.myID; myGroupID = settings.myGroupID;
    myFitTime = settings.fitTimeCrosscorr;
    % Execute delayed scatter script
    MW_delayedScatter
    
    % Renaming for later use
    concentrationCorrData = CorrData;
    concentrationassociatedFieldNames = associatedFieldNames;
    concentrationbadSchnitzes = badSchnitzes; 
    if exist('contourPlotData','var')
        concentrationContourPlotData = contourPlotData;
    end   
    
    % create correlation functions R(rate, growth)
    % ===
    
    % Set up appropriate field names for R(rate, mu)
    associatedFieldNames = {settings.timeFieldNameDerivative, settings.fluorDerivativeFieldName, settings.muFieldNameDerivative};
    % parse some settings from Excel file
    badSchnitzes = settings.badSchnitzes; alreadyRemovedInMatFile = settings.alreadyRemovedInMatFile;
    myID = settings.myID; myGroupID = settings.myGroupID;
    myFitTime = settings.fitTimeCrosscorr;
    % Execute delayed scatter script
    MW_delayedScatter
    
    % Renaming for later use
    rateCorrData = CorrData;
    rateassociatedFieldNames = associatedFieldNames;
    ratebadSchnitzes = badSchnitzes; 
    if exist('contourPlotData','var')
        rateContourPlotData = contourPlotData;
    end  
    
    % create autocorrelation functions R(Y, Y)
    % ===
    
    % Set up appropriate field names
    
    % growth
    associatedFieldNames = {settings.timeFieldName, settings.muFieldName, settings.muFieldName};
    MW_delayedScatter
    MW_autoCorr_corrtime_and_fitExponential
    % rename for later use
    growthautoCorrData = CorrData; growthautoFieldNames = associatedFieldNames; growthautoBadSchnitzes = badSchnitzes; 
    
    % concentration
    associatedFieldNames = {settings.timeFieldName, settings.fluorFieldName, settings.fluorFieldName};
    MW_delayedScatter
    MW_autoCorr_corrtime_and_fitExponential
    % rename for later use
    concentrationautoCorrData = CorrData; concentrationautoFieldNames = associatedFieldNames; concentrationautoBadSchnitzes = badSchnitzes; 
    
    % rate
    associatedFieldNames = {settings.timeFieldNameDerivative, settings.fluorDerivativeFieldName, settings.fluorDerivativeFieldName};
    MW_delayedScatter
    MW_autoCorr_corrtime_and_fitExponential        
    % rename for later use
    rateautoCorrData = CorrData; rateautoFieldNames = associatedFieldNames; rateautoBadSchnitzes = badSchnitzes; 
    
end    

%% If desired, save to appropriate database

% Additional options for preliminary analysis
% MW_summaryplotspreliminaryanalysis

% Additional options for full analysis
% ===

% config(..).m
% savingFluorDynamicsData

% For scripts directory see:
% - D:\Local_Software\Martijn_fluorDynamicsScripts

% For datasets see e.g. 
% - U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\CRPcAMP
% - \\storage01\data\AMOLF\groups\tans-group\Biophysics\2014-2017_Cellular-experiments\Fluor-validation\data








