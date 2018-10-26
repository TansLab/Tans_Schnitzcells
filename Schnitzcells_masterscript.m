
%% INSTRUCTIONS

% This file will allow you to perform all of the analyses that we generally
% perform to analyse growing microcolonies of cells. 
% 
% For each dataset, for each microcolony, an Excel file should be present
% that holds the parameter settings for the analysis. An example excel file
% should be present in the root directory. I suggest you copy this Excel
% file, rename it such that you can relate it to your dataset (current
% format Schnitzcells_Config_<experiment-date>_<microcolony-pos#>), and
% edit it to your needs.7
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
% (2)
% Look at configfile and DJK_initschnitz, w. respect to p.micronsPerPixel 
% and p.camera; if used outside our lab, this is not convenient..
% (3)
% fitTimeScatterHist in the config file seems to be an unused parameter.

%% Run this section to start the GUI interface
% (Press ctrl+enter when cursor is here.)

global ranFromGUI

if ~exist('runsections','var') || strcmp(runsections, 'none')

    runsections = 'none';
    GUIHandle = MW_GUI_schnitzcells
    ranFromGUI=1;

end

%% Parameters to set before starting.
% Set settings.analysisType and press CTRL+enter.
% The settings.analysisType parameter allows you to do a full analysis
% (standard), or preliminary analysis. The preliminary analysis does not
% track lineages, and allows you to get an impression of the whole colony
% behavior. (The preliminary analysis will create a schnitzcells file, but
% this file will not allow for a appropriate full analysis of your data.)
%
% ourSettings.analysisType     should be 'preliminary' or 'full' to
%                           automatically run the appropriate sections for
%                           a preliminary or full analysis.


% declaration of necessary global paramteters (for running from GUI)
global p ourSettings

% If not set yet, set flag that you don't run from GUI
if ~exist('ranFromGUI','var')
    ranFromGUI=0;
end

% In case you DO NOT want to run from GUI, the following code is
% automatically executed when runsections is not set by the GUI
if ~ranFromGUI
            
    % Parameters to set
    % ==
    
    % Full or prelinary analysis
    if ~isfield(ourSettings,'analysisType')
        % ourSettings.analysisType = 'preliminary'; % preliminary or full
        ourSettings.analysisType = 'full';
    end   

    % If not manually determined, set the runsections parameter appropriately
    if ~exist('runsections','var')
        if strcmp(ourSettings.analysisType, 'preliminary')
            runsections = 'allpreliminary'; % all is default
        elseif strcmp(ourSettings.analysisType, 'full')
            runsections = 'allfull';        
        end
    end
    
    % (re)load data at end of analysis from disk if set to 1.
    LOADDATASETATEND = 1;
       
end

% Parameters you should NOT change
% ===
% line where list of parameters starts in Excel file.
ourSettings.EXCELREADSTART = 14; 
% text editor that matlab will attempt to open to show you files.
ourSettings.MYTEXTEDITOR = 'notepad.exe'; 


%% Choose the Excel file with setting for your dataset.
% Press CTRL+enter.
% Run this section and a dialogue box will allow you to choose an Excel
% file. This section simply saves the path of the file you have chosen to
% ourSettings.configfilepath.

if any(strcmp(runsections,{'allpreliminary', 'allfull' 'loadfile'}))
    
    %select the configfile to run
    [ourSettings.myconfigfilename, ourSettings.mypathname] = uigetfile('*.xls;*.xlsx','Select a Schnitzcells configuration file');
    % full path to file
    ourSettings.configfilepath = [ourSettings.mypathname ourSettings.myconfigfilename];
        
    %old way
    %ourSettings.mypathname = 'F:\A_Tans1_step1_incoming_not_backed_up\2015-12-02\';
    %myfilename = 'Schnitzcells_Analysis_Config_2015_12_02_pos2col1.xlsx';

end

%% Read configuration settings from excel file.
% Press CTRL+enter.
% The function MW_readsettingsfromexcelfile reads in the Excel file, and
% puts all settings in a struct called "ourSettings".

if any(strcmp(runsections,{'allpreliminary', 'allfull', 'loadfile','reloadfile'}))
    % One can execute this section again to reload ourSettings 
    % Note that these are then not immediate parsed to the "p" struct.    
    
    % load ourSettings from excel file
    % ourSettings.configfilepath gives the path to MW_readsettingsfromexcelfile
    [ourSettings, alldata] = MW_readsettingsfromexcelfile(ourSettings)
    
    % Update field in GUI w. name of config file
    if ranFromGUI
        set(handles.thefilename,'String',[ourSettings.mypathname ourSettings.myconfigfilename]);
    end
end

%% Double check whether to continue
% Press CTRL+enter.
% This step is mainly to prevent accidentally running whole script 

if any(strcmp(runsections,{'allpreliminary', 'allfull'}))
        
    analysisTypes = {'PRELIMINARY','FULL'}

    myAnswer = questdlg(['Loaded ' [ourSettings.mypathname ourSettings.myconfigfilename] ', for ' ourSettings.analysisType ' analysis do you want to continue?'],'Confirmation required.','Yes','No','No');
    if strcmp(myAnswer, 'No') || strcmp(myAnswer, 'Cancel')
        error('Analysis aborted.');
    end

end

%% Now make a vector with parameter values that schnitzcells scripts can handle. (and some misc. other admin)
% Press CTRL+enter.
% Schnitzcells takes the struct "p" as input. This is basically a less
% extended version of the ourSettings struct. "p" is set here.
%
% Note that the settings for microscope and camera will affect the number
% of microns per pixels, which is an important setting in determining the
% length.
% The microns per pixel are hard coded in the DJK_initschnitz function. The
% values are also listed below for reference:
% hamamatsu camera, setup 1 (new camera, installed in 2014-11)
% p.micronsPerPixel = 0.0438;
% CoolSnap camera, setup 1 (old camera)
% p.micronsPerPixel = 0.04065;
% camera hamamatsu, setup 2 (always same camera)
% p.micronsPerPixel = 0.04312;


if any(strcmp(runsections,{'allpreliminary', 'allfull','createp'}))
    
    disp('Now creating ''p'' struct from ourSettings struct.');

    % Create the p vector which holds all parameters and is fed into, and also
    % returned by most functions of the schnitzcells analysis software.
    % Use "ourSettings" struct as a base.
    % DJK_initschnitz only checks for errors and adds a few parameters based on
    % the already given paramteres.
    % ===
    % TODO this can be done more elegantly (but note that existence of two
    % vectors, "ourSettings" and "p", allows user to update "ourSettings" vector 
    % intermediately.
    p = DJK_initschnitz(ourSettings.positionName,ourSettings.movieDate,'e.coli.amolf','rootDir',...
        ourSettings.rootDir, 'cropLeftTop',ourSettings.cropLeftTop, 'cropRightBottom',ourSettings.cropRightBottom,...
        'fluor1',ourSettings.fluor1,'fluor2',ourSettings.fluor2,'fluor3',ourSettings.fluor3,...
        'setup',ourSettings.setup,'softwarePackage',ourSettings.softwarePackage,'camera',ourSettings.camera);
        % NOTE: p.movieDate should only contain the date (e.g. 1998-10-10), 
        % whilst ourSettings.movieDate often also contains an explanatory 
        % suffix.(E.g. 1998-10-10_OldSchoolExperiment.) Therefor, 
        % p.movieDate = ourSettings.movieDate(1:10). 
        % This is a bit awkward, but since
        % there are already 50+ configuration files that use this
        % convention, I did not change it when I realized this was
        % inconvenient.
        % -MW, 2018
    
    % Set framerange according to analysis type
    if any(strcmp(ourSettings.analysisType,'preliminary')) % fast analysis
        ourSettings.currentFrameRange = ourSettings.frameRangePreliminary;
    elseif any(strcmp(ourSettings.analysisType,'full')) % full analysis
        ourSettings.currentFrameRange = ourSettings.frameRangeFull;
    end
    
end
%% Some more hard coded parameters
% Press CTRL+enter.
% Set an ouput directory (mainly for the preliminary analysis).

if any(strcmp(runsections,{'allpreliminary', 'allfull','createp'}))

    % Output directory for figures etc.    
    ourSettings.MYOUTPUTDIR = [p.dateDir 'outputSummary\'];
    
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
ourSettings.performCropping = strcmp(questdlg('Select whether you in general want to crop your dataset. (Answer will be saved to Excel config file.','Cropping','Yes','No','No'),'Yes');
% save preference to excel
performCroppingIndex = find(strcmp({alldata{:,1}},'performCropping'))+ourSettings.EXCELREADSTART-1; % find line w. cropRightBottom field.    
if isempty(performCroppingIndex)
   error('Could not fiend performCroppingIndex Excel config file field. Maybe parameter performCropping is not set?');
end
xlswrite([ourSettings.mypathname ourSettings.myconfigfilename],{num2str(ourSettings.performCropping)},['B' num2str(performCroppingIndex) ':B' num2str(performCroppingIndex) '']); % write value to it

if ourSettings.performCropping
    
    if strcmp(ourSettings.analysisType, 'preliminary')

        % Manually make sure image dir is correct
        % (This is done to accomodate cropping.)
        p.imageDir = [ourSettings.rootDir ourSettings.movieDate '\' ourSettings.positionName '\']
        
        % set new croparea and save crop area to excel file
        % ===
        
        % Ask to crop, and ask 
        %myAnswer = questdlg(['Start cropping? And save selection to Excel file (close it first)?'],'Confirmation required.','Save,use,crop','Crop using old','Abort!','Save,use,crop');        
        myAnswer = questdlg(['Do you want to select a crop area?'],'Crop area.','Yes','Crop using current ourSettings','Abort!','Yes');        
        % Only select new if desired
        if strcmp(myAnswer, 'Yes')            
            
            % Determine crop area
            [selectedLeftTop,selectedRightBottom] = MW_determinecroparea(p, ourSettings.frameRangePreliminary);
            
            ourSettings.cropLeftTop = selectedLeftTop;
            ourSettings.cropRightBottom = selectedRightBottom;

            ExcelcropLeftTopIndex = find(strcmp({alldata{:,1}},'cropLeftTop'))+ourSettings.EXCELREADSTART-1; % find line w. cropLeftTop field.
            xlswrite([ourSettings.mypathname ourSettings.myconfigfilename],{mat2str(ourSettings.cropLeftTop)},['B' num2str(ExcelcropLeftTopIndex) ':B' num2str(ExcelcropLeftTopIndex) '']); % write value to it
            ExcelcropRightBottomIndex = find(strcmp({alldata{:,1}},'cropRightBottom'))+ourSettings.EXCELREADSTART-1; % find line w. cropRightBottom field.    
            xlswrite([ourSettings.mypathname ourSettings.myconfigfilename],{mat2str(ourSettings.cropRightBottom)},['B' num2str(ExcelcropRightBottomIndex) ':B' num2str(ExcelcropRightBottomIndex) '']); % write value to it
        end

        % cropping itself
        if strcmp(myAnswer, 'Yes') | strcmp(myAnswer, 'Crop using current ourSettings')
            % Crop images
            DJK_cropImages_3colors(p, ourSettings.frameRangePreliminary, ourSettings.cropLeftTop, ...
                ourSettings.cropRightBottom, 'cropName', [ourSettings.positionName ourSettings.cropSuffix]);    
        end

        disp('Done cropping');
    
    elseif strcmp(ourSettings.analysisType, 'full')
    
        % Manually make sure image dir is correct
        % (This is done to accomodate cropping.)
        p.imageDir = [ourSettings.rootDir ourSettings.movieDate '\' ourSettings.positionName '\']

        DJK_cropImages_3colors(p, ourSettings.frameRangeFull, ourSettings.cropLeftTop, ...
            ourSettings.cropRightBottom, 'cropName', [ourSettings.positionName ourSettings.cropSuffix]);

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
if ourSettings.performCropping

    % =========================================================================
    % Load this setting later if you want to skip cropping
    % =========================================================================
    % (After loading ourSettings from excel file.)
    p = DJK_initschnitz([ourSettings.positionName ourSettings.cropSuffix],ourSettings.movieDate,'e.coli.amolf','rootDir',...
        ourSettings.rootDir, 'cropLeftTop',ourSettings.cropLeftTop, 'cropRightBottom',ourSettings.cropRightBottom,...
        'fluor1',ourSettings.fluor1,'fluor2',ourSettings.fluor2,'fluor3',ourSettings.fluor3,...
        'setup',ourSettings.setup,'softwarePackage',ourSettings.softwarePackage,'camera',ourSettings.camera)
    % =========================================================================

    disp('p (parameters) loaded');

end
end

%% Segmenation
% Press CTRL+enter.
% The function PN_segmoviephase_3colors will segment all images in the
% framerange provided by ourSettings.currentFrameRange.
% =========================================================================
% ADVANCED PARAMETER ourSettings that can be useful:
% p.useFullImage=1; % forces to use full image for segmentation
% p.overwrite=1; % overwrite existing files (to redo segmentation)
% p.customColonyCenter=[x,y]; % set center of colony to select ROI 
% =========================================================================

if any(strcmp(runsections,{'allpreliminary', 'allfull','segmentation'}))

    PN_segmoviephase_3colors(p,'segRange', ourSettings.currentFrameRange,'slices', ourSettings.slices,...
        'rangeFiltSize', ourSettings.rangeFiltSize,'maskMargin', ourSettings.maskMargin,'LoG_Smoothing',...
        ourSettings.LoG_Smoothing,'minCellArea', ourSettings.minCellArea,...
        'GaussianFilter', ourSettings.GaussianFilter,'minDepth', ourSettings.minDepth,'neckDepth', ourSettings.neckDepth);

    PN_copySegFiles(p,'segRange', ourSettings.currentFrameRange,'slices', ourSettings.slices,...
        'rangeFiltSize', ourSettings.rangeFiltSize,'maskMargin', ourSettings.maskMargin,'LoG_Smoothing',...
        ourSettings.LoG_Smoothing,'minCellArea', ourSettings.minCellArea,...
        'GaussianFilter', ourSettings.GaussianFilter,'minDepth', ourSettings.minDepth,'neckDepth', ourSettings.neckDepth);

    % Let user now section is done by sound
    mysound=load('gong'); sound(mysound.y);

    disp('Done with segmentation.');
   
end

%% To manually redo frames, execute this code manually or via GUI
% Press CTRL+enter.
% This section is generally not executed; only when a frame's segmentation
% has failed, you can execute this code that allows you to redo that frame
% with different parameters.

if any(strcmp(runsections,{'redosegforframe'}))       
    if ~ranFromGUI
        % EDIT PARAMETERS HERE IF MANUALLY REDOING
        TOREDOFRAME = 10;
        SLICESTEMPORARY = [2]; % default [1 2 3] instead of ourSettings.slices
        LOGSMOOTHINGTEMPORARY = 10; % default 2; instead of ourSettings.LoG_Smoothing
        MINDEPTHTEMP = 5;% default 5;
        MINCELLAREA = 250; % default 250      
        RANGEFILTSIZE          = 35; % default 35
        MASKMARGIN             = 5; % default 5
        GAUSSIANFILTER         = 5; % default 5 
        NECKDEPTH              = 2; % default 2
    else    
        % ourSettings from GUI
        TOREDOFRAME            = ourSettings.TOREDOFRAME;
        SLICESTEMPORARY        = ourSettings.SLICESTEMPORARY;
        LOGSMOOTHINGTEMPORARY  = ourSettings.LOGSMOOTHINGTEMPORARY;
        MINDEPTHTEMP           = ourSettings.MINDEPTHTEMP;
        MINCELLAREA            = ourSettings.MINCELLAREA;
        RANGEFILTSIZE          = ourSettings.RANGEFILTSIZE;
        MASKMARGIN             = ourSettings.MASKMARGIN;
        GAUSSIANFILTER         = ourSettings.GAUSSIANFILTER;
        NECKDEPTH              = ourSettings.NECKDEPTH;
    end
    
    theOriginalFrameRange = ourSettings.currentFrameRange; ourSettings.currentFrameRange = TOREDOFRAME;
    p.overwrite=1;    
    
    % a set problemCells field functions as flag, so should 
    if isfield(p, 'problemCells'), p=rmfield(p, 'problemCells'), end
    
    PN_segmoviephase_3colors(p,'segRange', ourSettings.currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', RANGEFILTSIZE,'maskMargin', MASKMARGIN,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', GAUSSIANFILTER,'minDepth', MINDEPTHTEMP,'neckDepth', NECKDEPTH);

    PN_copySegFiles(p,'segRange', ourSettings.currentFrameRange,'slices', SLICESTEMPORARY,...
        'rangeFiltSize', RANGEFILTSIZE,'maskMargin', MASKMARGIN,'LoG_Smoothing',...
        LOGSMOOTHINGTEMPORARY,'minCellArea', MINCELLAREA,...
        'GaussianFilter', GAUSSIANFILTER,'minDepth', MINDEPTHTEMP,'neckDepth', NECKDEPTH);    
   
    ourSettings.currentFrameRange = theOriginalFrameRange;
    p.overwrite=0;
end

disp('Section done');

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
    if strcmp(ourSettings.analysisType, 'preliminary')
        ourSettings.assistedYesNo=0;
    elseif strcmp(ourSettings.analysisType, 'full')
        ourSettings.assistedYesNo=1;
    end

    % run function which allows user to manually check segmentation
    PN_manualcheckseg(p,'manualRange',ourSettings.currentFrameRange,'overwrite',0,'assistedCorrection',ourSettings.assistedYesNo); % ADVANCED PARAMETERS in fn arguments

end
%% Perform quick analysis of segmentation
% Press CTRL+enter.
% If you're doing a preliminary analysis, you can run DJK_analyzeSeg to get
% an impression of the growth behavior of the colony.

if any(strcmp(runsections,{'allpreliminary', 'allfull','quickanalysis'}))
    
    DJK_analyzeSeg(p,'manualRange',ourSettings.currentFrameRange,'onscreen',1,'DJK_saveDir',ourSettings.MYOUTPUTDIR);
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
    DJK_trackcomplete(p,'trackRange',ourSettings.currentFrameRange,'trackMethod','singleCell');    
    
    disp('Done tracking');
elseif any(strcmp(runsections,{'allfull','trackandmanualcorrections'}))% full
    p.overwrite=0; % ADVANCED SETTING
    p.showAll = 1; % show all segmented frames to user

    % Default ourSettings for identifying problem cells, when not given config
    % file
    if ~isfield(ourSettings,'pixelsMoveDef')
        ourSettings.pixelsMoveDef=15;
    end
    if ~isfield(ourSettings,'pixelsLenDef')
        ourSettings.pixelsLenDef=[-4 13];
    end
    if ~isfield(ourSettings,'pixelsAreaDiv')
        ourSettings.pixelsAreaDiv=[-4 13];
    end
    
    % Alternative trackers one can try
    % This is quite useful if errors occur. manualrange can be edited to
    % two frames only. Set p.overwrite to 1. Then re-run this whole 
    % section to update general tracking file. (p.overwrite is reset
    % automatically above).
    % ===
    % default tracker
    %DJK_tracker_djk(p,'manualRange', ourSettings.currentFrameRange)
    % slow but more robust:
    %NW_tracker_centroid_vs_area(p,'manualRange', [1:244]); 
    % fast but fails often:
    %MW_tracker(p,'manualRange', [1:244]); 
    % If all else fails:
    % edit MW_helperforlinkingframes
        
    DJK_tracker_djk(p,'manualRange', ourSettings.currentFrameRange); % default tracker           
    
    % Find problem cells
    [problems, theOutputFilePath] = DJK_analyzeTracking(p,'manualRange', ourSettings.currentFrameRange, ...
        'pixelsMoveDef', ourSettings.pixelsMoveDef, 'pixelsLenDef', ourSettings.pixelsLenDef,...
        'pixelsAreaDiv', ourSettings.pixelsAreaDiv);
    % open output file in external editor (not necessary, but convenient)
    eval(['!' ourSettings.MYTEXTEDITOR ' ' theOutputFilePath ' &']);
    
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
    PN_manualcheckseg(p,'manualRange',ourSettings.currentFrameRange,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted

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
% Set ourSettings.specialtracker to run a desired alternative tracker (either 
% MW or NW). If this field is not set, default tracker is run.

if any(strcmp(runsections,{'customtrackersoncustomrange'}))
   
    if ~isfield(p, 'overwrite')
        p.overwrite=0; % ADVANCED SETTING
    end
    if ~isfield(p, 'showAll');
        p.showAll = 1; % show all segmented frames to user
    end
    
    % call desired tracker
    if ~isfield(ourSettings, 'specialtracker') 
        DJK_tracker_djk(p,'manualRange', ourSettings.retrackFrameRange); % default tracker
    elseif strcmp(ourSettings.specialtracker, 'MW')
        MW_tracker(p,'manualRange', ourSettings.retrackFrameRange); 
    elseif strcmp(ourSettings.specialtracker, 'NW')
        NW_tracker_centroid_vs_area(p,'manualRange', ourSettings.retrackFrameRange);
    end
    
    % Now, if retracked range was not full range, the overall tracking
    % should be redone. This function is already executed in the trackers,
    % but should be re-executed here in this case.
    if ~isequal(ourSettings.retrackFrameRange, ourSettings.currentFrameRange)
        disp('CREATING SCHNITZCELLS FROM FULL DESIRED RANGE ***');
        MW_calculateSchnitzPropertiesWithoutTracking(p,'manualRange', ourSettings.currentFrameRange);
    end
    
end


%% (Optional) Check again after alternative tracker
% <TODO>
% If in previous section alternative tracker was called, the tracking also
% needs to be checked for suspcious cells. Do this with the earlier
% section. This section is under construction.

if any(strcmp(runsections,{'checkaftercustom'}))
    
    % Default ourSettings for identifying problem cells, when not given config
    % file
    if ~isfield(ourSettings,'pixelsMoveDef')
        ourSettings.pixelsMoveDef=15;
    end
    if ~isfield(ourSettings,'pixelsLenDef')
        ourSettings.pixelsLenDef=[-4 13];
    end
    if ~isfield(ourSettings,'pixelsAreaDiv')
        ourSettings.pixelsAreaDiv=[-4 13];
    end
    if ~isfield(p,'overwrite')
        p.overwrite=0;
    end

    %disp('Option not working yet.. Use (re)track (..) instead.');
    
    
    % Find problem cells    
    [problems, theOutputFilePath] = DJK_analyzeTracking(p,'manualRange', ourSettings.currentFrameRange, ...
        'pixelsMoveDef', ourSettings.pixelsMoveDef, 'pixelsLenDef', ourSettings.pixelsLenDef,...
        'pixelsAreaDiv', ourSettings.pixelsAreaDiv);
    % open output file in external editor (not necessary, but convenient)
    eval(['!' ourSettings.MYTEXTEDITOR ' ' theOutputFilePath ' &']);
    
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
    PN_manualcheckseg(p,'manualRange',ourSettings.currentFrameRange,'override',p.overwrite,'assistedCorrection',0); % assisted correction of because problem cells highlighted

    
    
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
    DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1,'manualRange',ourSettings.currentFrameRange);

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

%%
if any(strcmp(runsections,{'allfull','makemovieraw'}))

    % user info
    disp('Making movie, manual options like p.hideFig and p.showNr can be found in ''help MW_makeMovieRaw''');
    
    if ~isfield(p,'showNr')
        p.showNr=2;
    end
    if ~isfield(p,'hideFig')
        p.hideFig=1;
    end
    
    % Make movie
    MW_makeMovieRaw(p,ourSettings);

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
    NW_initializeFluorData(p,'manualRange', ourSettings.currentFrameRange);
    % Load PSF (color-independent)
    load(ourSettings.fluorPointSpreadFunctionPath, 'PSF');

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
        load(ourSettings.fluorCorrectionImagePaths{colorIdx}, 'flatfield', 'shading', 'replace');

        
        % Finding shifts.
        disp('Looking for shifts');
        if ~isfield(p,'mothermachine')
            optimalShift = DJK_getFluorShift_anycolor(p,'manualRange', ourSettings.currentFrameRange,'fluorcolor',currentFluor,'maxShift',MAXSHIFT);
                % BUG HERE? TODO: the same shift is applied to all images. I
                % wonder whether this is what is intended..
            disp('Correcting');
        else
            disp('INFO: mothermachine option is activated, so not looking for optimal fluor shift. Set optimalShift = [0,0].')
            optimalShift = [0,0];            
        end
        
        % Correct images (shading, background).
        DJK_correctFluorImage_anycolor(p, flatfield, shading, replace,'manualRange', ourSettings.currentFrameRange,  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF),'fluorcolor',currentFluor,'minimalMode',0);

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
    
    % In case of full analysis, also calculate skeleton length
    if strcmp(ourSettings.analysisType, 'full')
        % Also calculate skeleton lengths
        NDL_addToSchnitzes_skeletonLengthMW(p); % saves schnitzcells param
        % And make sure DJK_addToSchnitzes_mu also uses length_skeleton
        p.lengthFields = {'rp_length' 'length_fitCoef3b' 'length_fitNew' 'length_skeleton'};         
    end

    %8 Add correct mu
    DJK_addToSchnitzes_mu(p, 'frameSizes', ourSettings.muWindow );% p.lengthFields is set above

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
            % Here the field "phase2" is added, which is similar to phase, 
            % but contains also "predicted" values for 
            % for fields where phase was not known. Predictions are based
            % on another parameter (length, in this case).
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
        % also for new rate:      
        schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
            ['d' X '5_divAreaPx'],... % dX5
            ['phase2_at_d' X '5_time']);
        % and for concentration
        schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
            [X '6_mean'],... % X6mean
            ['phase2_at' X '']); % etc..

        for muWindowSize = ourSettings.muWindow
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
        for muWindowSize = ourSettings.muWindow

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
        
        % Now also add a field that identifies frames at which
        % concentration and production rate of fluor have been determiend.
        p.schnitzFilePath = [p.tracksDir p.movieName '-Schnitz.mat'];
        schnitzcells = MW_addToSchnitzes_indices_atXdX(p.schnitzFilePath,X);
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
    [fitTime, fitMu] = DJK_analyzeMu(p, schnitzcells, 'onScreen', 1,'fitTime',[0 10000],'DJK_saveDir',ourSettings.MYOUTPUTDIR);
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
        [fitFluorMean(colorIdx), fitFluorVariance(colorIdx)] = DJK_plot_avColonyOverTime(p, schnitzcells, fluorFieldName, 'fitTime', fitTime, 'onScreen', 1,'DJK_saveDir',ourSettings.MYOUTPUTDIR);
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
    xlswrite([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{'Identifier'},['A1:A1'])
    xlswrite([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{[p.movieDate '_' p.movieName]},['A' lineToWriteTo ':' 'A' lineToWriteTo])
    
    for i = 1:numel(summaryParameters)
        letterToWriteTo = ['' i+64+1]; % +1 since start at B
        %disp(['writing ' [letterToWriteTo lineToWriteTo]]);
        xlswrite([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],{summaryParametersNames{i}},[letterToWriteTo '1:' letterToWriteTo '1'])
        xlswrite([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.xls'],summaryParameters(i),[letterToWriteTo lineToWriteTo ':' letterToWriteTo lineToWriteTo])
    end

    % Save output to .mat file (update the matfile)
    if exist([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'file') == 2
        load([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata');
    end
    thedata(posNumber).summaryParameters = summaryParameters;
    thedata(posNumber).ourSettings = ourSettings;
    thedata(posNumber).p = p; % "ourSettings" and "p" are a bit redundant for historic reasons.
    save([ourSettings.MYOUTPUTDIR 'summaryParametersPreliminary.mat'],'thedata','summaryParametersNames');
    
disp('Done making summary preliminary analysis.');
disp('Check out MW_summaryplotspreliminaryanalysis.m for more plots when all positions done');
    
end


%% DONE / LOADING DATA
% Press CTRL+enter.
% Load data if desired.

if any(strcmp(runsections,{'allpreliminary', 'allfull', 'makeoutputfull','rerunfullanalysis'}))

    % We are done with ourSettings up schnitzcells structure! All data is saved.
    % You can proceed with analyzing the data..

    % To load after full analysis is done:
    if exist('LOADDATASETATEND','var'), if LOADDATASETATEND

        % Load this dataset using
        % =========================================================================
        p = DJK_initschnitz([ourSettings.positionName ourSettings.cropSuffix], ourSettings.movieDate,'e.coli.AMOLF','rootDir',...
         ourSettings.rootDir, 'cropLeftTop', ourSettings.cropLeftTop, 'cropRightBottom', ourSettings.cropRightBottom,...
             'fluor1',ourSettings.fluor1,...
             'fluor2',ourSettings.fluor2,...
             'fluor3',ourSettings.fluor3,...
             'setup',ourSettings.setup,...
             'softwarePackage',ourSettings.softwarePackage,...
             'camera',ourSettings.camera);

        [p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);
        % =========================================================================

    end, end

    disp('Loading complete.');

end

%% Adding the production rate field
if any(strcmp(runsections,{'allfull', 'makeoutputfull', 'rerunfullanalysis'})) % full

% This is backwards compatibility to add the field dX5_divAreaPx and
% dX5_divAreaPx_cycCor.

% First gather fluor color letters
allFluorColors = {p.fluor1 p.fluor2 p.fluor3};
fluorColorsUsed = {};
for flIdx = 1:3
    if ~strcmp(allFluorColors{flIdx},'none')
        fluorColorsUsed{end+1} = upper(allFluorColors{flIdx});
    end
end
% then determine additional fields if necessary
saveFlag=0;
for flIdx = 1:numel(fluorColorsUsed)
    X = fluorColorsUsed{flIdx};

    if ~isfield(schnitzcells, ['d' X '5_divAreaPx'])
        schnitzcells = MW_fluorRate_anycolor(schnitzcells,X,'length_fitNew');

        % Note that the new rate is compatible (only) with X5 phase
        schnitzcells = DJK_addToSchnitzes_cycleCor(schnitzcells,...
            ['d' X '5_divAreaPx'],... % dX5
            ['phase2_at_d' X '5_time']);

        saveFlag=1;
    end
end

% Save the result
if saveFlag
    NW_saveSchnitzcells(p,schnitzcells); 
end

end
%% More analyses for full analysis
% Press CTRL+enter.
% This section provides a load of analysis of the schnitzcells structure.
% It calculates noise values and cross-correlations of important
% parameters. It is convenient to save this analysis to a database 
% with all your experiments. The script savingFluorDynamicsData provides
% you with this opportunity (see section "Database" below).

if any(strcmp(runsections,{'allfull', 'makeoutputfull', 'rerunfullanalysis'})) % full
    
    %% Seting up main output parameter
    output = struct;
    
    % ourSettings up some more parameters
    p.myID = ourSettings.myID
    ourSettings.myOutputFolder = [ourSettings.mypathname  '\' p.movieDate  '_' p.movieName '_' p.myID  '\'];
    
    p.NW_saveDir = [ourSettings.myOutputFolder 'misc\'];  % To send additional output to
    p.DJK_saveDir = [ourSettings.myOutputFolder 'misc\']; % To send additional output to
    
    % Location of .mat file containing schnitzcells struct
    ourSettings.myDataFile = [ourSettings.mypathname '\' p.movieName  '\data\' p.movieName '-Schnitz.mat'];
   
    % Load datafile
    load(ourSettings.myDataFile);

    % Some more parameter renaming
    myFile = ourSettings.myDataFile; EXPORTFOLDER = ourSettings.myOutputFolder; FITTIME = ourSettings.fitTimeMu;    
    
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
    ourSettings.theLetters='';
    for fluorIdx = 1:3
        % check whether this fluor field was set by user in config file
        if isfield(ourSettings,['fluor' num2str(fluorIdx)])
        if ~strcmp(upper(ourSettings.(['fluor' num2str(fluorIdx)])),'NONE')
            % determine fluor code letter
            theLetter = upper(ourSettings.(['fluor' num2str(fluorIdx)]));
            ourSettings.theLetters = [ourSettings.theLetters theLetter];
            % set fields
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName = ...
                strrep(ourSettings.timeFieldName,'X',theLetter);
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName = ...
                strrep(ourSettings.fluorFieldName,'X',theLetter);
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldName = ...
                strrep(ourSettings.muFieldName,'X',theLetter)
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative = ...
                strrep(ourSettings.timeFieldNameDerivative,'X',theLetter);
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName = ...
                strrep(ourSettings.fluorDerivativeFieldName,'X',theLetter);
            ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldNameDerivative = ...
                strrep(ourSettings.muFieldNameDerivative,'X',theLetter);
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
        PLOTSCATTER=ourSettings.PLOTSCATTER; % activate to plot all scatter plots
        
        % Set up appropriate field names for R(concentration,mu)
        % TODO MW: todo: make this X_time etc, and do a strrep for X to applicable color
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldName};
        
        % obtain some ourSettings from Excel file
        badSchnitzes = ourSettings.badSchnitzes; alreadyRemovedInMatFile = ourSettings.alreadyRemovedInMatFile;
        addSlowOnes = ourSettings.addSlowOnes;
        myID = ourSettings.myID; myGroupID = ourSettings.myGroupID;
        myFitTime = ourSettings.fitTimeCrosscorr;
        myOutputFolder = ourSettings.myOutputFolder;
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
        PLOTSCATTER=ourSettings.PLOTSCATTER; % activate to plot all scatter plots
        
        % Set up appropriate field names for R(rate, mu)
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldNameDerivative};
        % parse some ourSettings from Excel file
        badSchnitzes = ourSettings.badSchnitzes; alreadyRemovedInMatFile = ourSettings.alreadyRemovedInMatFile;
        addSlowOnes = ourSettings.addSlowOnes;
        myID = ourSettings.myID; myGroupID = ourSettings.myGroupID;
        myFitTime = ourSettings.fitTimeCrosscorr;
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
        addSlowOnes=0; % the third field is not mu any more
        
        % Set up appropriate field names

        % growth
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).muFieldName};
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
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName};
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
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName};
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
    
        %% R(prod, E): cross-correlations within same fluor between production and concentation
        % This might be particularly usefull for detecting feedback
        PLOTSCATTER=0;
        addSlowOnes=0; % the third field is not mu any more
        
        % The concentration field has more datapoint than the rate
        % field, due the nature of how the derivative is calculated.
        
        % Therefor, gather the applicable concentration field values based
        % on the time vectors of the dX field.
        someSuffix   = ['_at_d' ourSettings.theLetters(fluorIdx)];
        if ~isfield(schnitzcells, [ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName '_at_d' ourSettings.theLetters(fluorIdx)])
            %%
            disp('Preparing additional concentration field');
            schnitzcellsPath    = ourSettings.myDataFile;
            timeFieldSelection  = ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative;
            timeFieldSource     = ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName;
            sourceField         = ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName;            
            schnitzcells = MW_addToSchnitzes_fieldYatTime(schnitzcellsPath,timeFieldSelection, timeFieldSource, sourceField, someSuffix);
        end
        
        associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, [ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName someSuffix], ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName};
        % execute analysis scripts
        MW_delayedScatter
        % rename for later use
        output.rateConcentrationCC{fluorIdx} = CorrData; 
        output.rateConcentrationCCFieldNames{fluorIdx} = associatedFieldNames; 
        output.rateConcentrationCCBadschnitzes{fluorIdx} = badSchnitzes;                 
        
        
        %% Create cross-correlation functions for fluorescent colors
        disp('Making crosscorrs for color vs. color');
        
        if numel(activeFluorIdx)>1
        for fluorIdx2 = fluorIdx+1:max(activeFluorIdx)            
            
            dualCounter=dualCounter+1; 
            
            %% concentration fluor N vs fluor M
            associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx2).fluorFieldName};
            % execute analysis scripts
            MW_delayedScatter
            MW_autoCorr_corrtime_and_fitExponential
            % rename for later use
            output.concentrationDualCrossCorrData{dualCounter}          = CorrData; 
            output.concentrationDualCrossCorrFieldNames{dualCounter}    = associatedFieldNames; 
            output.concentrationDualCrossCorrBadSchnitzes{dualCounter}  = badSchnitzes; 

            %% rate fluor N vs fluor M
            associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx2).fluorDerivativeFieldName};
            % execute analysis scripts
            MW_delayedScatter
            MW_autoCorr_corrtime_and_fitExponential        
            % rename for later use
            output.rateDualCrossCorrData{dualCounter}       = CorrData; 
            output.rateDualCrossFieldNames{dualCounter}     = associatedFieldNames; 
            output.rateDualCrossBadSchnitzes{dualCounter}   = badSchnitzes; 
            
            %% rate fluor N vs concentration fluor M
            % Note that we can look at 
            % Concentration_N -> production_N (make sense to observe feedback system)
            % Production_N -> Concentration_N
            % Concentration_N -> production_M (makes sense for ribosomes)
            % Concentration_M -> production_N
            % Production_M -> Concentration_N (doesn't seem applicable yet)
            % Production_N -> Concentration_M
            
            % TODO: Add code below to analyze the N,M correlations (see
            % above for correlations of between the same protein in
            % production and concentration).
            
            %{
            associatedFieldNames = {ourSettings.fieldNamesWithFluorLetter(fluorIdx).timeFieldNameDerivative, ourSettings.fieldNamesWithFluorLetter(fluorIdx).fluorDerivativeFieldName, ourSettings.fieldNamesWithFluorLetter(fluorIdx2).fluorFieldName};                    
            
            % execute analysis scripts
            MW_delayedScatter
            MW_autoCorr_corrtime_and_fitExponential        
            % rename for later use
            output.rateDualCrossCorrData{dualCounter}       = CorrData; 
            output.rateDualCrossFieldNames{dualCounter}     = associatedFieldNames; 
            output.rateDualCrossBadSchnitzes{dualCounter}   = badSchnitzes; 
            %}
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
    winopen(ourSettings.myOutputFolder);
    
end    

%% Save p and ourSettings struct to file in output dir if desired
if any(strcmp(runsections,{'allfull', 'makeoutputfull','rerunfullanalysis'})) % full
    
   if ~exist('DONTSAVEFULLANALYSISATEND','var')
    
       outputFilename = [p.dateDir 'outputandsettings_v2_' p.movieDate '_' p.movieName '.mat'];
       save(outputFilename, 'p', 'ourSettings', 'schnitzcells', 's_rm', 'output');

       disp(['p (parameter), output, schnitzcells and ourSettings structs were saved to: ' outputFilename]);
       
   end
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








