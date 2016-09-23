
%% INSTRUCTIONS

% This file will allow you to perform all of the analyses that we generally
% perform to analyse growing microcolonies of cells. 
% 
% For each dataset, for each microcolony, an Excel file should be present
% that holds the parameter settings for the analysis. An example excel file
% should be present in the root directory. I suggest you copy this Excel
% file, rename it such that you can relate it to your dataset (current
% format Schnitzcells_Config_<experiment-date>_<microcolony-pos#>), and
% edit it to your needs.
%
% If you have edited your Excel file appropriately, continue with the
% analysis below by executing the script section by section. The parameter
% "runsections", which is set automatically, decides which sections are
% executed.
% In fact, if you execute the script without any settings, a GUI will open
% that allows you to automatically run the script section by section. The
% matlab file associated with this GUI is:
% >> MW_GUI_schnitzcells 
% The main output of the analysis (aside from plots) is the data file which
% contains the "schnitzcells" struct. Once generated, it will be located 
% in pos1\data\pos1-Schnitz.mat. It contains the lineage structure with all
% the separate individual cells (i.e. birth-to-division, aka one 
% "schnitz").
% 
% 
% A more extensive description of what the matlab scripts do can be found
% in the comments on each section.
%
% MW, 2016/03
% 
% ===
%
% This software was originally developed in the Elowitz group, see:
% Young, J. W., Locke, J. C. W., Altinok, A., Rosenfeld, N., Bacarian, T., 
% Swain, P. S., … Elowitz, M. B. (2012). Measuring single-cell gene 
% expression dynamics in bacteria using fluorescence time-lapse microscopy. 
% Nature Protocols, 7(1), 80–8. doi:10.1038/nprot.2011.432
%
% It has been extensively edited and extended in the tans biophysics lab.
% It is currently stored at (private repository):
% https://bitbucket.org/microscopeguerrillas/schnitzcells_tans
% To request access contact the Tans Biophysics lab,
% http://www.amolf.nl/research/biophysics/
% http://tansgroup.amolf.nl/
% Martijn "MW" Wehrens, wehrens@amolf.nl; PI: Sander Tans, tans@amolf.nl.


%% TODOs

% (1)
% Note sure whether preliminary analysis should save as much to the
% schnitzcells .mat file as it does now. Maybe it is smarter to separate
% the info that analyis saves to a separate file. (This might have been
% implemented better in the excell files with command lists that were
% previously used..)


%% Run this section to start the GUI interface
% (Press ctrl+enter when cursor is here.)

if ~exist('runsections','var') || strcmp(runsections, 'none')

    runsections = 'none';
    MW_GUI_schnitzcells

end

%% Parameters to set before starting.
% Set settings.analysisType and press CTRL+enter.
% The settings.analysisType parameter allows you to do a full analysis
% (standard), or preliminary analysis. The preliminary analysis does not
% track lineages, and allows you to get an impression of the whole colony
% behavior. (The preliminary analysis will create a schnitzcells file, but
% this file will not allow for a appropriate full analysis of your data.)
%
% settings.analysisType     should be 'preliminary' or 'full' to
%                           automatically run the appropriate sections for
%                           a preliminary or full analysis.


% declaration of necessary global paramteters (for running from GUI)
global p settings

% In case you DO NOT want to run from GUI, the following code is
% automatically executed when runsections is not set by the GUI
if ~exist('runsections', 'var')
            
    % Parameters to set
    % ==
    
    % Full or prelinary analysis
    if ~isfield(settings,'analysisType')
        % settings.analysisType = 'preliminary'; % preliminary or full
        settings.analysisType = 'full';
    end   

    % Set the runsections parameter appropriately
    if strcmp(settings.analysisType, 'preliminary')
        runsections = 'allpreliminary'; % all is default
    elseif strcmp(settings.analysisType, 'full')
        runsections = 'allfull';        
    end
    
    % (re)load data at end of analysis from disk if set to 1.
    LOADDATASETATEND = 1;
    
    % Signal to script whether has been run from GUI
    ranFromGUI = 0;
else
    ranFromGUI = 1;
end

% Parameters you should NOT change
% ===
% line where list of parameters starts in Excel file.
settings.EXCELREADSTART = 14; 
% text editor that matlab will attempt to open to show you files.
settings.MYTEXTEDITOR = 'notepad.exe'; 


%% Choose the Excel file with settings for your dataset.
% Press CTRL+enter.
% Run this section and a dialogue box will allow you to choose an Excel
% file. This section simply saves the path of the file you have chosen to
% settings.configfilepath.

if any(strcmp(runsections,{'allpreliminary', 'allfull' 'loadfile'}))
    
    %select the configfile to run
    [settings.myconfigfilename, settings.mypathname] = uigetfile('*.xls;*.xlsx','Select a Schnitzcells configuration file');
    % full path to file
    settings.configfilepath = [settings.mypathname settings.myconfigfilename];
        
    %old way
    %settings.mypathname = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-12-02\';
    %myfilename = 'Schnitzcells_Analysis_Config_2015_12_02_pos2col1.xlsx';

end

%% Read configuration settings from excel file.
% Press CTRL+enter.
% The function MW_readsettingsfromexcelfile reads in the Excel file, and
% puts all settings in a struct called "settings".

if any(strcmp(runsections,{'allpreliminary', 'allfull', 'loadfile','reloadfile'}))
    % One can execute this section again to reload settings 
    % Note that these are then not immediate parsed to the "p" struct.    
    
    % load settings from excel file
    % settings.configfilepath gives the path to MW_readsettingsfromexcelfile
    [settings, alldata] = MW_readsettingsfromexcelfile(settings)
    
    % Update field in GUI w. name of config file
    if ranFromGUI
        set(handles.thefilename,'String',[settings.mypathname settings.myconfigfilename]);
    end
end

%% Double check whether to continue
% Press CTRL+enter.
% This step is mainly to prevent accidentally running whole script 

if any(strcmp(runsections,{'allpreliminary', 'allfull'}))
        
    analysisTypes = {'PRELIMINARY','FULL'}

    myAnswer = questdlg(['Loaded ' [settings.mypathname settings.myconfigfilename] ', for ' settings.analysisType ' analysis do you want to continue?'],'Confirmation required.','Yes','No','No');
    if strcmp(myAnswer, 'No') || strcmp(myAnswer, 'Cancel')
        error('Analysis aborted.');
    end

end

%% Now make a vector with parameter values that schnitzcells scripts can handle. (and some misc. other admin)
% Press CTRL+enter.
% Schnitzcells takes the struct "p" as input. This is basically a less
% extended version of the settings struct. "p" is set here.

if any(strcmp(runsections,{'allpreliminary', 'allfull','createp'}))
    
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
        'setup',settings.setup,'softwarePackage',settings.softwarePackage,'camera',settings.camera);

    % Set framerange according to analysis type
    if any(strcmp(settings.analysisType,'preliminary')) % fast analysis
        settings.currentFrameRange = settings.frameRangePreliminary;
    elseif any(strcmp(settings.analysisType,'full')) % full analysis
        settings.currentFrameRange = settings.frameRangeFull;
    end
    
end
%% Some more hard coded parameters
% Press CTRL+enter.
% Set an ouput directory (mainly for the preliminary analysis).

if any(strcmp(runsections,{'allpreliminary', 'allfull','createp'}))

    % Output directory for figures etc.    
    settings.MYOUTPUTDIR = [p.dateDir 'outputSummary\'];
    
end

%% Crop images
% Press CTRL+enter.
% This section allows you to crop the images. If you don't want to crop the
% images, you need to manually copy images from the appropriate data 
% directory to a subdirecotry called ./images/.

% Crop images 
% =========================================================================

% Crop images
% Puts cropped images in new directory, which gets the suffix 
if any(strcmp(runsections,{'allpreliminary','cropimages'}))

% Ask user whether cropping is desired.
settings.performCropping = strcmp(questdlg('Select whether you in general want to crop your dataset. (Answer will be saved to Excel config file.','Cropping','Yes','No','No'),'Yes');
% save preference to excel
performCroppingIndex = find(strcmp({alldata{:,1}},'performCropping'))+settings.EXCELREADSTART-1; % find line w. cropRightBottom field.    
if isempty(performCroppingIndex)
   error('Could not fiend performCroppingIndex Excel config file field. Maybe parameter performCropping is not set?');
end
xlswrite([settings.mypathname settings.myconfigfilename],{num2str(settings.performCropping)},['B' num2str(performCroppingIndex) ':B' num2str(performCroppingIndex) '']); % write value to it

if settings.performCropping
    
    if strcmp(settings.analysisType, 'preliminary')

        % Manually make sure image dir is correct
        % (This is done to accomodate cropping.)
        p.imageDir = [settings.rootDir settings.movieDate '\' settings.positionName '\']
        
        % set new croparea and save crop area to excel file
        % ===
        
        % Ask to crop, and ask 
        %myAnswer = questdlg(['Start cropping? And save selection to Excel file (close it first)?'],'Confirmation required.','Save,use,crop','Crop using old','Abort!','Save,use,crop');        
        myAnswer = questdlg(['Do you want to select a crop area?'],'Crop area.','Yes','Crop using current settings','Abort!','Yes');        
        % Only select new if desired
        if strcmp(myAnswer, 'Yes')            
            
            % Determine crop area
            [selectedLeftTop,selectedRightBottom] = MW_determinecroparea(p, settings.frameRangePreliminary);
            
            settings.cropLeftTop = selectedLeftTop;
            settings.cropRightBottom = selectedRightBottom;

            ExcelcropLeftTopIndex = find(strcmp({alldata{:,1}},'cropLeftTop'))+settings.EXCELREADSTART-1; % find line w. cropLeftTop field.
            xlswrite([settings.mypathname settings.myconfigfilename],{mat2str(settings.cropLeftTop)},['B' num2str(ExcelcropLeftTopIndex) ':B' num2str(ExcelcropLeftTopIndex) '']); % write value to it
            ExcelcropRightBottomIndex = find(strcmp({alldata{:,1}},'cropRightBottom'))+settings.EXCELREADSTART-1; % find line w. cropRightBottom field.    
            xlswrite([settings.mypathname settings.myconfigfilename],{mat2str(settings.cropRightBottom)},['B' num2str(ExcelcropRightBottomIndex) ':B' num2str(ExcelcropRightBottomIndex) '']); % write value to it
        end

        % cropping itself
        if strcmp(myAnswer, 'Yes') | strcmp(myAnswer, 'Crop using current settings')
            % Crop images
            DJK_cropImages_3colors(p, settings.frameRangePreliminary, settings.cropLeftTop, ...
                settings.cropRightBottom, 'cropName', [settings.positionName settings.cropSuffix]);    
        end

        disp('Done cropping');
    
    elseif strcmp(settings.analysisType, 'full')
    
        % Manually make sure image dir is correct
        % (This is done to accomodate cropping.)
        p.imageDir = [settings.rootDir settings.movieDate '\' settings.positionName '\']

        DJK_cropImages_3colors(p, settings.frameRangeFull, settings.cropLeftTop, ...
            settings.cropRightBottom, 'cropName', [settings.positionName settings.cropSuffix]);

        disp('Done cropping');
        
    else
        
        disp('Did nothing.');
        
    end
    
end  
end


%% Update p accordingly.
% Press CTRL+enter.
% If you have cropped the images the p struct needs to be updated with the
% new location of the images. 

if any(strcmp(runsections,{'allpreliminary', 'allfull','cropimages','loadpforcropped'}))
if settings.performCropping

    % =========================================================================
    % Load this setting later if you want to skip cropping
    % =========================================================================
    % (After loading settings from excel file.)
    p = DJK_initschnitz([settings.positionName settings.cropSuffix],settings.movieDate,'e.coli.amolf','rootDir',...
        settings.rootDir, 'cropLeftTop',settings.cropLeftTop, 'cropRightBottom',settings.cropRightBottom,...
        'fluor1',settings.fluor1,'fluor2',settings.fluor2,'fluor3',settings.fluor3,...
        'setup',settings.setup,'softwarePackage',settings.softwarePackage,'camera',settings.camera)
    % =========================================================================

    disp('p (parameters) loaded');

end
end

%% Segmenation
% Press CTRL+enter.
% The function PN_segmoviephase_3colors will segment all images in the
% framerange provided by settings.currentFrameRange.
% =========================================================================
% ADVANCED PARAMETER SETTINGS that can be useful:
% p.useFullImage=1; % forces to use full image for segmentation
% p.overwrite=1; % overwrite existing files (to redo segmentation)
% p.customColonyCenter=[x,y]; % set center of colony to select ROI 
% =========================================================================

if any(strcmp(runsections,{'allpreliminary', 'allfull','segmentation'}))

    PN_segmoviephase_3colors(p,'segRange', settings.currentFrameRange,'slices', settings.slices,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

    PN_copySegFiles(p,'segRange', settings.currentFrameRange,'slices', settings.slices,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        settings.LoG_Smoothing,'minCellArea', settings.minCellArea,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', settings.minDepth,'neckDepth', settings.neckDepth);

    % Let user now section is done by sound
    mysound=load('gong'); sound(mysound.y);

    disp('Done with segmentation.');
   
end

%% To manually redo 1 frame, execute this code manually or via GUI
% Press CTRL+enter.
% This section is generally not executed; only when a frame's segmentation
% has failed, you can execute this code that allows you to redo that frame
% with different parameters.

if any(strcmp(runsections,{'redosegforframe'}))   
    if ~ranFromGUI
        % EDIT PARAMETERS HERE IF MANUALLY REDOING
        TOREDOFRAME = 10;
        SLICESTEMPORARY = [2]; % default [1 2 3] instead of settings.slices
        LOGSMOOTHINGTEMPORARY = 10; % default 2; instead of settings.LoG_Smoothing
        MINDEPTHTEMP = 5;% default 5;
        MINCELLAREA = 250; % default 250        
    else    
        % Settings from GUI
        TOREDOFRAME            = settings.TOREDOFRAME;
        SLICESTEMPORARY        = settings.SLICESTEMPORARY;
        LOGSMOOTHINGTEMPORARY  = settings.LOGSMOOTHINGTEMPORARY;
        MINDEPTHTEMP           = settings.MINDEPTHTEMP;
        MINCELLAREA            = settings.MINCELLAREA;
    end
    
    theOriginalFrameRange = settings.currentFrameRange; settings.currentFrameRange = TOREDOFRAME;
    p.overwrite=1;    
    
    % a set problemCells field functions as flag, so should 
    if isfield(p, 'problemCells'), p=rmfield(p, 'problemCells'), end
    
    PN_segmoviephase_3colors(p,'segRange', settings.currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', MINDEPTHTEMP,'neckDepth', settings.neckDepth);

    PN_copySegFiles(p,'segRange', settings.currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', settings.rangeFiltSize,'maskMargin', settings.maskMargin,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', settings.GaussianFilter,'minDepth', MINDEPTHTEMP,'neckDepth', settings.neckDepth);    
   
    settings.currentFrameRange = theOriginalFrameRange;
    p.overwrite=0;
end



%% Start the manual checking by user
% Press CTRL+enter.
% The segmentation performed by the Matlab script often contains errors.
% There is an extensive script to interact with the segmentation files from
% each frame, in order to correct errors in the segmentation. The is done
% by the function PN_manualcheckseg, which will provide a list of commands
% that allow you to navigate and edit through the segmentation files.
% The segmentation files can be found in the ./segmentation directory, and
% contain the raw segmentation (LNsub) plus - after running
% PN_manualcheckseg - the corrected segmentation that is used for the
% analysis in the matrix Lc.

if any(strcmp(runsections,{'allpreliminary', 'allfull','manualchecksegfull'}))

    % choose option based type analysis
    if strcmp(settings.analysisType, 'preliminary')
        settings.assistedYesNo=0;
    elseif strcmp(settings.analysisType, 'full')
        settings.assistedYesNo=1;
    end

    % run function which allows user to manually check segmentation
    PN_manualcheckseg(p,'manualRange',settings.currentFrameRange,'overwrite',0,'assistedCorrection',settings.assistedYesNo); % ADVANCED PARAMETERS in fn arguments

end
%% Perform quick analysis of segmentation
% Press CTRL+enter.
% If you're doing a preliminary analysis, you can run DJK_analyzeSeg to get
% an impression of the growth behavior of the colony.

if any(strcmp(runsections,{'allpreliminary', 'allfull','quickanalysis'}))
    
    DJK_analyzeSeg(p,'manualRange',settings.currentFrameRange,'onscreen',1,'DJK_saveDir',settings.MYOUTPUTDIR);
    %close(gcf);
    
    disp('Quick analysis done.');
    
end

%% Perform tracking, check tracking, make problem movie and correct
% Press CTRL+enter. To reiterate: press CTRL+enter.
% Iterate this section until no problems are listed any more in:
% \analysis\tracking\manualRangeX_Y\posZcrop-tracking.txt
% This analysis will track cells from frame to frame (using the default
% function DJK_tracker_djk, see next section for alternatives), thus
% identifying lineages. However, often there are mistakes in this procedure. 
% These can often be corrected by correcting the segmentation. Therefor,
% this section also tries to identify suspicious events in the lineage, and
% suspicous cells are listed in aformentioned text file. Consequtively, the
% segmentation correction can be run again. When p.problemCells is set it
% will highlight the suspicious cells listed there. One can press "s" to
% give frame or schnitzcells cell numbers. When all frames are manually
% corrected, one can run this section again to see whether all problems
% have disappeared. (Note that suspicious cells not always turn out to be
% wrongly tracked.)
% If you have identified schnitzes (indivual cells), or lineages, that show
% issues (e.g. artefacts that cannot be removed), the cells can be marked
% as "bad schnitzes" in the Excel file. Lineages, or lineages containing
% these cells, will later be excluded from the analysis.
% Here, for the first time the schnitzcells file is generated, which
% contains all the output of the analyses. This file can be fuond in
% posX\data\posX-Schnitz.mat

if any(strcmp(runsections,{'allpreliminary','trackpreliminary'})) % fast
    DJK_trackcomplete(p,'trackRange',settings.currentFrameRange,'trackMethod','singleCell');    
    
    disp('Done tracking');
elseif any(strcmp(runsections,{'allfull','trackandmanualcorrections'}))% full
    p.overwrite=0; % ADVANCED SETTING
    p.showAll = 1; % show all segmented frames to user

    % Default settings for identifying problem cells, when not given config
    % file
    if ~isfield(settings,'pixelsMoveDef')
        settings.pixelsMoveDef=15;
    end
    if ~isfield(settings,'pixelsLenDef')
        settings.pixelsLenDef=[-4 13];
    end
    
    % Alternative trackers one can try
    % This is quite useful if errors occur. manualrange can be edited to
    % two frames only. Set p.overwrite to 1. Then re-run this whole 
    % section to update general tracking file. (p.overwrite is reset
    % automatically above).
    % ===
    % default tracker
    %DJK_tracker_djk(p,'manualRange', settings.currentFrameRange)
    % slow but more robust:
    %NW_tracker_centroid_vs_area(p,'manualRange', [1:244]); 
    % fast but fails often:
    %MW_tracker(p,'manualRange', [1:244]); 
    % If all else fails:
    % edit MW_helperforlinkingframes
        
    DJK_tracker_djk(p,'manualRange', settings.currentFrameRange); % default tracker           
    
    % Find problem cells
    [problems, theOutputFilePath] = DJK_analyzeTracking(p,'manualRange', settings.currentFrameRange, 'pixelsMoveDef', settings.pixelsMoveDef, 'pixelsLenDef', settings.pixelsLenDef);
    % open output file in external editor (not necessary, but convenient)
    eval(['!' settings.MYTEXTEDITOR ' ' theOutputFilePath ' &']);
    
    p.problemCells = problems;

    % Create lookup table for frame label <-> schnitz label
    p.slookup = MW_makeslookup(p);

    % Let user know we're done
    mysound=load('gong'); sound(mysound.y);

    % Manual check again 
    % Since "p" now contains the lookup table, and problemcells is defined, it
    % will highlight the problemcells.
    % Note that when one performs manual corrections to a frame, the mapping is
    % not correct any more.
    PN_manualcheckseg(p,'manualRange',settings.currentFrameRange,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted

    disp('Done (full) tracking');
    
end



%% (Optional) Some alternative trackers
% Press CTRL+enter. Re-iterate previous section when done.
% The code contains two alternative trackers, the defualt tracker is
% DJK_tracker_djk, but two alternatives are MW_tracker (fast, works by
% overlap, high failure rate if cells move a load), and
% NW_tracker_centroid_vs_area (slow, but often offers good tracking where 
% DJK_tracker_djk has failed).
% This section will only be run if customtrackersoncustomrange is set.
% Set settings.specialtracker to run a desired alternative tracker (either 
% MW or NW). If this field is not set, default tracker is run.

if any(strcmp(runsections,{'customtrackersoncustomrange'}))
   
    if ~isfield(p, 'overwrite')
        p.overwrite=0; % ADVANCED SETTING
    end
    if ~isfield(p, 'showAll');
        p.showAll = 1; % show all segmented frames to user
    end
    
    % call desired tracker
    if ~isfield(settings, 'specialtracker') 
        DJK_tracker_djk(p,'manualRange', settings.retrackFrameRange); % default tracker
    elseif strcmp(settings.specialtracker, 'MW')
        MW_tracker(p,'manualRange', settings.retrackFrameRange); 
    elseif strcmp(settings.specialtracker, 'NW')
        NW_tracker_centroid_vs_area(p,'manualRange', settings.retrackFrameRange);
    end
    
    % Now, if retracked range was not full range, the overall tracking
    % should be redone. This function is already executed in the trackers,
    % but should be re-executed here in this case.
    if ~isequal(settings.retrackFrameRange, settings.currentFrameRange)
        disp('CREATING SCHNITZCELLS FROM FULL DESIRED RANGE ***');
        MW_calculateSchnitzPropertiesWithoutTracking(p,'manualRange', settings.currentFrameRange);
    end
    
end


%% (Optional) Check again after alternative tracker
% <TODO>
% If in previous section alternative tracker was called, the tracking also
% needs to be checked for suspcious cells. Do this with the earlier
% section. This section is under construction.

if any(strcmp(runsections,{'checkaftercustom'}))
    
    % Default settings for identifying problem cells, when not given config
    % file
    if ~isfield(settings,'pixelsMoveDef')
        settings.pixelsMoveDef=15;
    end
    if ~isfield(settings,'pixelsLenDef')
        settings.pixelsLenDef=[-4 13];
    end
    
    %disp('Option not working yet.. Use (re)track (..) instead.');
    
    
    % Find problem cells    
    [problems, theOutputFilePath] = DJK_analyzeTracking(p,'manualRange', settings.currentFrameRange, 'pixelsMoveDef', 15, 'pixelsLenDef', [-4 13]);
    % open output file in external editor (not necessary, but convenient)
    eval(['!' settings.MYTEXTEDITOR ' ' theOutputFilePath ' &']);
    
    p.problemCells = problems;

    % Create lookup table for frame label <-> schnitz label
    p.slookup = MW_makeslookup(p);

    % Let user know we're done
    mysound=load('gong'); sound(mysound.y);

    % Manual check again 
    % Since "p" now contains the lookup table, and problemcells is defined, it
    % will highlight the problemcells.
    % Note that when one performs manual corrections to a frame, the mapping is
    % not correct any more.
    PN_manualcheckseg(p,'manualRange',settings.currentFrameRange,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted

    
    
end

%% Once above section has been iterated to satisfaction make sure to remove problemCells field from p
% Press CTRL+enter after having reiterated earlier steps to satisfaction.
% (Because movies later also work with this.)

if any(strcmp(runsections,{'allfull','makemovie'}))
    
    % remove problemCells parameter from p if it exists (because otherwise
    % the function that makes the movie will assume you're only interested
    % in those frames.)
    if isfield(p,'problemCells');
        p =rmfield(p,'problemCells');
    end
end

%% Make movie that allows you to double check the segmentation & tracking
% Press CTRL+enter.
% This sections creates a movie of the growing microcolony, with extra
% (color) information about lineages. Check it to see whether all mistakes
% in identifying lineages have been removed. The movie is exported to the
% analysis directory, subdirectory ./movies, but this folder is opened
% automatically when script is done.

if any(strcmp(runsections,{'allfull','makemovie'}))

    % User info
    if isfield(p, 'problemCells')
        disp('p.problemCells field exists, will make movie of those frames only.');
    end

    % Make movie
    DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1,'manualRange',settings.currentFrameRange);

    %{
    % Other movie making options
    DJK_makeMovie(p, 'tree', 'schAll', 'stabilize', 1,'problemCells',problems);
    DJK_makeMovie (p, 'tree', 'cellno', 'stabilize', 1,'problemCells',problems);
    p=rmfield(p,'problemCells')
    DJK_makeMovie (p, 'tree', 'cellno', 'stabilize', 1);
    %DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1);
    %}
    
    % Let user know done, open output folder.
    mysound=load('gong'); sound(mysound.y);    
    disp(['Done making analysis. Opening movie folder at ' p.analysisDir 'movies\..'])
    winopen([p.analysisDir 'movies\']);
end

%% Create backup of segtrack
% Press CTRL+enter.
% Creates a backup of the tracking.

if any(strcmp(runsections,{'allfull','createbackup'}))
    NW_Backup_SegTrack(p);
end

%% Fluor part I: Start of correcting fluor
% Press CTRL+enter.
% Some initializing functions are executed here, but the most important
% function is DJK_correctFluorImage_anycolor, which corrects fluor images
% with shading, background, etc.

if any(strcmp(runsections,{'allpreliminary', 'allfull','correctionsandanalysis'}))
    MAXSHIFT = 10;

    % Initialization, applies to all colors
    % ===
    NW_initializeFluorData(p,'manualRange', settings.currentFrameRange);
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
        optimalShift = DJK_getFluorShift_anycolor(p,'manualRange', settings.currentFrameRange,'fluorcolor',currentFluor,'maxShift',MAXSHIFT);
        disp('Correcting');
        % Correct images (shading, background).
        DJK_correctFluorImage_anycolor(p, flatfield, shading, replace,'manualRange', settings.currentFrameRange,  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF),'fluorcolor',currentFluor,'minimalMode',0);

    end

    % Add fluor information to schnitz struct (schnitzcells struct is saved).
    disp('Adding fluor info to schnitzcells structure.')
    DJK_compileSchnitzImproved_3colors(p,'quickMode',0);

    % Let user hear we're done
    mysound=load('gong'); sound(mysound.y);

    disp('Done correcting fluor part I.');
end

%% Lengths: Now calculate lengths and add them to schnitzcells struct
% Press CTRL+enter.
% Functions executed here calculate the length and growth rate, of all of 
% the cells in the microcolony. DJK_addToSchnitzes_length calculates
% length, DJK_addToSchnitzes_mu calculates growth rate.
% Some of the fluorescence analysis also needs the fields generated here

if any(strcmp(runsections,{'allpreliminary', 'allfull','correctionsandanalysis'}))

    % ADVANCED PARAMETERS
    %p.schnitzNum='all';
    %p.onscreen=1;

    DJK_addToSchnitzes_length(p);
    
    % In case of full analysis
    if strcmp(settings.analysisType, 'full')
        % Also calculate skeleton lengths
        NDL_addToSchnitzes_skeletonLengthMW(p);
        % And make sure DJK_addToSchnitzes_mu also uses length_skeleton
        p.lengthFields = {'rp_length' 'length_fitCoef3b' 'length_fitNew' 'length_skeleton'};         
    end

    %8 Add correct mu
    DJK_addToSchnitzes_mu(p, 'frameSizes', settings.muWindow); % p.lengthFields is set above

    % Let user know we're done
    mysound=load('gong'); sound(mysound.y);
    disp('Done calculating lengths.');

end
%% Fluor part II: Continue with fluor again now
% Press CTRL+enter.
% Add values of fluor concentration, fluor production to the schnitzcells
% struct.

if any(strcmp(runsections,{'allpreliminary', 'allfull','correctionsandanalysis'}))

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

    disp('Done with all preparing analyses. Now you can analyse data.');
end

%% For preliminary analysis, create some data output
% Press CTRL+enter.
% Setting up the schnitzcells struct is done now. 
% If you were doing a preliminary analysis, this
% section will provide some analysis. 
%  ===

if any(strcmp(runsections,{'allpreliminary','analysispreliminary'})) 

    % Obtain schnitzcells
    [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

    % Fit mu
    [fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'onScreen', 1,'fitTime',[0 10000],'DJK_saveDir',settings.MYOUTPUTDIR);
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
        [fitFluorMean(colorIdx), fitFluorVariance(colorIdx)] = DJK_plot_avColonyOverTime(p, schnitzcells, fluorFieldName, 'fitTime', fitTime, 'onScreen', 1,'DJK_saveDir',settings.MYOUTPUTDIR);
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
    xlswrite([settings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{'Identifier'},['A1:A1'])
    xlswrite([settings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{[p.movieDate '_' p.movieName]},['A' lineToWriteTo ':' 'A' lineToWriteTo])
    
    for i = 1:numel(summaryParameters)
        letterToWriteTo = ['' i+64+1]; % +1 since start at B
        %disp(['writing ' [letterToWriteTo lineToWriteTo]]);
        xlswrite([settings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{summaryParametersNames{i}},[letterToWriteTo '1:' letterToWriteTo '1'])
        xlswrite([settings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],summaryParameters(i),[letterToWriteTo lineToWriteTo ':' letterToWriteTo lineToWriteTo])
    end

    % Save output to .mat file (update the matfile)
    if exist([settings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'file') == 2
        load([settings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata');
    end
    thedata(posNumber).summaryParameters = summaryParameters;
    thedata(posNumber).settings = settings;
    thedata(posNumber).p = p; % "settings" and "p" are a bit redundant for historic reasons.
    save([settings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata','summaryParametersNames');
    
disp('Done making summary preliminary analysis.');
disp('Check out MW_summaryplotspreliminaryanalysis.m for more plots when all positions done');
    
end


%% DONE / LOADING DATA
% Press CTRL+enter.
% Load data if desired.

if any(strcmp(runsections,{'allpreliminary', 'allfull', 'makeoutputfull','rerunfullanalysis'}))

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

    disp('Loading complete.');

end
%% More analyses for full analysis
% Press CTRL+enter.
% This section provides a load of analysis of the schnitzcells structure.
% It calculates noise values and cross-correlations of important
% parameters. It is convenient to save this analysis to a database 
% with all your experiments. The script savingFluorDynamicsData provides
% you with this opportunity (see section "Database" below).



if any(strcmp(runsections,{'allfull', 'makeoutputfull', 'rerunfullanalysis'})) % full
    
    % Seting up main output parameter
    output = struct;
    
    % Settings up some more parameters
    p.myID = settings.myID
    settings.myOutputFolder = [settings.mypathname  '\' p.movieDate  '_' p.movieName '_' p.myID  '\'];
    
    p.NW_saveDir = [settings.myOutputFolder 'misc\'];  % To send additional output to
    p.DJK_saveDir = [settings.myOutputFolder 'misc\']; % To send additional output to
    
    % Location of .mat file containing schnitzcells struct
    settings.myDataFile = [settings.mypathname '\' p.movieName  '\data\' p.movieName '-Schnitz.mat'];
   
    % Load datafile
    load(settings.myDataFile);

    % Some more parameter renaming
    myFile = settings.myDataFile; EXPORTFOLDER = settings.myOutputFolder; FITTIME = settings.fitTimeMu;    
    
    % Set up export directory
    if ~exist(EXPORTFOLDER,'dir'), mkdir(EXPORTFOLDER); end
    
    % Some plots of cell growth
    % ===
    % execute script
    MW_growthplots
    
    %% Create correlation parameter names
    % ===
    activeFluorIdx = [];
    % loop over different fluor slots
    settings.theLetters='';
    for fluorIdx = 1:3
        % check whether this fluor field was set by user in config file
        if isfield(settings,['fluor' num2str(fluorIdx)])
        if ~strcmp(upper(settings.(['fluor' num2str(fluorIdx)])),'NONE')
            % determine fluor code letter
            theLetter = upper(settings.(['fluor' num2str(fluorIdx)]));
            settings.theLetters = [settings.theLetters theLetter];
            % set fields
            settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName = ...
                strrep(settings.timeFieldName,'X',theLetter);
            settings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName = ...
                strrep(settings.fluorFieldName,'X',theLetter);
            settings.fieldNamesWithFluorLetter(fluorIdx).muFieldName = ...
                strrep(settings.muFieldName,'X',theLetter)
            settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative = ...
                strrep(settings.timeFieldNameDerivative,'X',theLetter);
            settings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName = ...
                strrep(settings.fluorDerivativeFieldName,'X',theLetter);
            settings.fieldNamesWithFluorLetter(fluorIdx).muFieldNameDerivative = ...
                strrep(settings.muFieldNameDerivative,'X',theLetter);
            % mark as active
            activeFluorIdx(end+1)=fluorIdx;
        end  
        end
    end
    
    %% Create all plots for all fluors
    dualCounter=0;
    for fluorIdx = activeFluorIdx
        

        %% create correlation functions R(concentration, growth)
        % ===
        PLOTSCATTER=settings.PLOTSCATTER; % activate to plot all scatter plots
        
        % Set up appropriate field names for R(concentration,mu)
        % TODO MW: todo: make this X_time etc, and do a strrep for X to applicable color
        associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).muFieldName};
        % obtain some settings from Excel file
        badSchnitzes = settings.badSchnitzes; alreadyRemovedInMatFile = settings.alreadyRemovedInMatFile;
        myID = settings.myID; myGroupID = settings.myGroupID;
        myFitTime = settings.fitTimeCrosscorr;
        myOutputFolder = settings.myOutputFolder;
        % Execute delayed scatter script
        MW_delayedScatter


        %% Renaming for later use
        % output.branchavg{fieldname} contains average data for branches
        output.concentrationBranch_groups{fluorIdx} = branch_groups; 
        output.concentrationBranch_groupsControl{fluorIdx} = branch_groupsControl;
        
        output.concentrationCorrData{fluorIdx} = CorrData; 
        output.concentrationCorrDataControl{fluorIdx} = CorrDataControl;
        
        output.concentrationassociatedFieldNames{fluorIdx} = associatedFieldNames;
        output.concentrationbadSchnitzes{fluorIdx} = badSchnitzes; 
        if exist('contourPlotData','var')
            output.concentrationContourPlotData{fluorIdx} = contourPlotData;
        end    

        %% create correlation functions R(rate, growth)
        % ===
        PLOTSCATTER=settings.PLOTSCATTER; % activate to plot all scatter plots
        
        % Set up appropriate field names for R(rate, mu)
        associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, settings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).muFieldNameDerivative};
        % parse some settings from Excel file
        badSchnitzes = settings.badSchnitzes; alreadyRemovedInMatFile = settings.alreadyRemovedInMatFile;
        myID = settings.myID; myGroupID = settings.myGroupID;
        myFitTime = settings.fitTimeCrosscorr;
        % Execute delayed scatter script
        MW_delayedScatter

        % Renaming for later use
        output.rateCorrData{fluorIdx} = CorrData; 
        output.rateCorrDataControl{fluorIdx} = CorrDataControl;
        
        output.rateassociatedFieldNames{fluorIdx} = associatedFieldNames;
        output.ratebadSchnitzes{fluorIdx} = badSchnitzes;     
        if exist('contourPlotData','var')
            output.rateContourPlotData{fluorIdx} = contourPlotData;
        end  

        %% create autocorrelation functions R(Y, Y)
        % ===
        PLOTSCATTER=0; % activate to plot all scatter plots
        
        % Set up appropriate field names

        % growth
        associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).muFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).muFieldName};
        % execute analysis scripts
        MW_delayedScatter
        MW_autoCorr_corrtime_and_fitExponential
        % rename for later use
        output.growthautoCorrData = CorrData; 
        output.growthautoFieldNames = associatedFieldNames; 
        output.growthautoBadSchnitzes = badSchnitzes;     
        % note that noise parameters are based on all branch data lumped
        % together, it would be better to take raw schnitzcells data as input.
        output.growthNoise = theNoise; 
        output.growthMean = theMean; 
        output.growthStd = theStd;

        %% concentration
        associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName};
        % execute analysis scripts
        MW_delayedScatter
        MW_autoCorr_corrtime_and_fitExponential
        % rename for later use
        output.concentrationautoCorrData{fluorIdx} = CorrData; 
        output.concentrationautoFieldNames{fluorIdx} = associatedFieldNames; 
        output.concentrationautoBadSchnitzes{fluorIdx} = badSchnitzes; 
        % note that noise parameters are based on all branch data lumped
        % together, it would be better to take raw schnitzcells data as input.
        output.concentrationNoise{fluorIdx} = theNoise; 
        output.concentrationMean{fluorIdx} = theMean; 
        output.concentrationStd{fluorIdx} = theStd;

        %% rate
        associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, settings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName};
        % execute analysis scripts
        MW_delayedScatter
        MW_autoCorr_corrtime_and_fitExponential        
        % rename for later use
        output.rateautoCorrData{fluorIdx} = CorrData; 
        output.rateautoFieldNames{fluorIdx} = associatedFieldNames; 
        output.rateautoBadSchnitzes{fluorIdx} = badSchnitzes; 
        % note that noise parameters are based on all branch data lumped
        % together, it would be better to take raw schnitzcells data as input.
        output.rateNoise{fluorIdx} = theNoise; 
        output.rateMean{fluorIdx} = theMean; 
        output.rateStd{fluorIdx} = theStd;
    
        %% Create cross-correlation functions for fluorescent colors
        disp('Making crosscorrs for color vs. color');
        
        if numel(activeFluorIdx)>1
        for fluorIdx2 = fluorIdx+1:max(activeFluorIdx)            
            
            dualCounter=dualCounter+1; 
            
            %% concentration fluor N vs fluor M
            associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, settings.fieldNamesWithFluorLetter(fluorIdx2).fluorFieldName};
            % execute analysis scripts
            MW_delayedScatter
            MW_autoCorr_corrtime_and_fitExponential
            % rename for later use
            output.concentrationDualCrossCorrData{dualCounter}          = CorrData; 
            output.concentrationDualCrossCorrFieldNames{dualCounter}    = associatedFieldNames; 
            output.concentrationDualCrossCorrBadSchnitzes{dualCounter}  = badSchnitzes; 

            %% rate fluor N vs fluor M
            associatedFieldNames = {settings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, settings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, settings.fieldNamesWithFluorLetter(fluorIdx2).fluorDerivativeFieldName};
            % execute analysis scripts
            MW_delayedScatter
            MW_autoCorr_corrtime_and_fitExponential        
            % rename for later use
            output.rateDualCrossCorrData{dualCounter}       = CorrData; 
            output.rateDualCrossFieldNames{dualCounter}     = associatedFieldNames; 
            output.rateDualCrossBadSchnitzes{dualCounter}   = badSchnitzes; 
            
        end
        end
    end
    
    %% Clear savedirs in case reran
    if isfield(p,'NW_saveDir')
        p=rmfield(p,'NW_saveDir');
    end
    if isfield(p,'DJK_saveDir')
        p=rmfield(p,'DJK_saveDir');
    end
    
    %% Let user know done, open output folder.
    disp('Done making analysis. Opening output folder.')
    winopen(settings.myOutputFolder);
    
end    

%% Save p and settings struct to file in output dir if desired
if any(strcmp(runsections,{'allfull', 'makeoutputfull','rerunfullanalysis'})) % full
   outputFilename = [p.dateDir 'outputandsettings_' p.movieName '.mat'];
   save(outputFilename, 'p', 'settings', 'schnitzcells', 's_rm', 'output');
   
   disp(['p (parameter), output, schnitzcells and settings structs were saved to: ' outputFilename]);
end


%% Database. If desired, save to appropriate database

% Additional options for preliminary analysis
% MW_summaryplotspreliminaryanalysis

% Additional options for full analysis
% ===

% CONFIGFILE =  'config(..).m'
% savingFluorDynamicsData

% For scripts directory see:
% - D:\Local_Software\Martijn_fluorDynamicsScripts

% For datasets see e.g. 
% - U:\ZZ_EXPERIMENTAL_DATA\A_Step5_Data_per_project\CRPcAMP
% - \\storage01\data\AMOLF\groups\tans-group\Biophysics\2014-2017_Cellular-experiments\Fluor-validation\data








