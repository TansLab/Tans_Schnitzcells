

% FUNCTION IS DEPRECATED
%
% Use option p1.useFullImage = 1 with PN_segmoviephase_3colors instead.
% 
% This function will be removed.
% - MW 2015/04

















function p = NW_nocut_segmoviephase_3colors_diffStrains(p,varargin)
%
% OPTIONAL ARGUMENTS:
% 'medium'    : 'rich'  . Special preparation of images. default: 'normal'
% (18.1.12: obsolete, nonfunctional! NW)
% 'method'    : 'brightfield'. Pretreatment of images in French style
%               (brightfield and 32 layers. 1 image per layer). default:
%               'phasecontrast'

% numimagesperlayer : which of the repeated images of the same layer are to
%               be taken and averaged. Default: all. Possible alternative: [1]
% quickMode: 1: Load directly averaged images from previous treatment.
%            Default: 0

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking. Use inputParser in
% recent matlab releases.
%-------------------------------------------------------------------------------

disp('WARNING: This function is deprecated, use p1.useFullImage = 1 with PN_segmoviephase_3colors instead.');

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

%%%%%%%%%%%default values
%STEP A : finds a global mask and crop the image
if ~existfield(p,'rangeFiltSize')                         %typical area for dectection of interesting features of the image
    p.rangeFiltSize = 35;
end
if ~existfield(p,'maskMargin')                            %additional margin of the mask : enlarge if cells missing on the edge
    p.maskMargin = 20;
end
%STEP B : find edges
if ~existfield(p,'LoG_Smoothing')                         %smoothing amplitude of the edge detection filter
    p.LoG_Smoothing = 2;
end
if ~existfield(p,'minCellArea')                           %minimum cell area (objects smaller than that will be erased)
    p.minCellArea = 250;
end
%STEP C : prepare seeds for watershedding
if ~existfield(p,'GaussianFilter')                        %smoothing of the original image to find local minima within cells
    p.GaussianFilter = 5;
end
if ~existfield(p,'minDepth')                              %minimum accepted depth for a local minimum
    p.minDepth = 5;
end
%STEP E: treatment of long cells
if ~existfield(p,'neckDepth')                             %minimum neck width to cut a too long cell
    p.neckDepth = 2;
end
%saving images
if ~existfield(p,'saveSteps')                             %indicate if you want to save intermediate images
    p.saveSteps = true;
end
if ~existfield(p,'PN_saveDir') & p.saveSteps              %subfolder of p.segmentationDir where image treatment steps are saved
    p.PN_saveDir = ['param' '_Marg' num2str(p.maskMargin) '_LoG' num2str(p.LoG_Smoothing) '_Area' num2str(p.minCellArea) '_Depth' num2str(p.minDepth) '_Neck' num2str(p.neckDepth) filesep];
    [message errmsg] = sprintf(['Image saving folder automatically generated: ' p.PN_saveDir]);
    disp(message);
end
%growth medium (nb: special problems in rich medium: holes + fringed edges)
if ~existfield(p,'medium')                                %special treatment in case of rich medium
    p.medium='normal';  
end
if strcmp(p.medium,'rich')==0 && strcmp(p.medium,'normal')==0
    disp('unknown kind of medium. Set to ''normal'' ');
    p.medium='normal';
end
%acquisition technique: 'phasecontrast' (conventional) or 'brightfield' (French
%style) possible
if ~existfield(p,'method')                                %standard technique
    p.method='phasecontrast';  
end
if strcmp(p.method,'phasecontrast')==0 && strcmp(p.method,'brightfield')==0
    disp('unknown kind of acquisition method. Set to ''phasecontrast'' ');
    p.method='phasecontrast';
end
%saveDir for averaged BrightField Images
if ~existfield(p,'NW_brightfieldsaveDir') %subfolder of p.imageDir
    p.NW_brightfieldsaveDir = ['averagedimages'];
end
%saveDir for correlated BrightField Images
if ~existfield(p,'NW_correlatedsaveDir') %subfolder of p.imageDir
    p.NW_correlatedsaveDir = ['correlatedimages'];
end
%quickMode
if ~existfield(p,'quickMode') %extract already averaged images
    p.quickMode = 0;
end

%--------------------------------------------------------------------------
% Checking or creation of directories
%--------------------------------------------------------------------------
% make sure every directory field has a trailing filesep
if (p.PN_saveDir(end) ~= filesep)
    p.PN_saveDir = [p.PN_saveDir filesep];
end
% if directory doesn't exist, create it
if exist(p.PN_saveDir)~=7
    [status,msg,id] = mkdir([p.segmentationDir p.PN_saveDir]);
    if status == 0
        disp(['Warning: unable to mkdir ' p.PN_saveDir ' : ' msg]);
    end
end
% Als make sure \png\ directory exists
png_dir = [p.segmentationDir 'png' filesep];
if exist(png_dir)~=7
    [status,msg,id] = mkdir([png_dir]);
    if status == 0
        disp(['Warning: unable to mkdir ' png_dir ' : ' msg]);
        return;
    end
end
if (strcmp(p.method,'brightfield')==1)
    % create extra directory for average images (brightfield only) (NW 23.1.12)
    if (p.NW_brightfieldsaveDir(end) ~= filesep)
        p.NW_brightfieldsaveDir = [p.NW_brightfieldsaveDir filesep];
    end
    if exist(p.NW_brightfieldsaveDir)~=7
        [status,msg,id] = mkdir([p.imageDir p.NW_brightfieldsaveDir]);
        if status == 0
            disp(['Warning: unable to mkdir ' p.NW_brightfieldsaveDir ' : ' msg]);
        end
    end
    % create extra directory for correlated images (brightfield only) (NW 23.1.12)
    if (p.NW_correlatedsaveDir(end) ~= filesep)
        p.NW_correlatedsaveDir = [p.NW_correlatedsaveDir filesep];
    end
    if exist(p.NW_correlatedsaveDir)~=7
        [status,msg,id] = mkdir([p.imageDir p.NW_correlatedsaveDir]);
        if status == 0
            disp(['Warning: unable to mkdir ' p.NW_correlatedsaveDir ' : ' msg]);
        end
    end
end

%-------------------------------------------------------------------------------
% Figure out which images exist
%--------------------------------------------------------------------------

% figure out the number of image layers (slices) and the number of images (in case of brightfield)
% per position capture

%--------------------------------------------------------------------------
% PHASECONTRAST
if (strcmp(p.method,'phasecontrast')==1)
    Dphase1   = dir([p.imageDir, '*-p-1-*.tif']);
    DphaseAll = dir([p.imageDir, '*-p-*.tif']);
    if ~isfield(p,'numphaseslices')
        if isempty(DphaseAll)
            error(['Can''t find any images in directory ' p.imageDir]);
        else
            p.numphaseslices = int8(ceil(length(DphaseAll)/length(Dphase1)));
            disp(['You appear to have ' num2str(p.numphaseslices) ' phase image per frame.']);
        end
    end

    if ~isfield(p,'slices')
        p.slices = [ 1:p.numphaseslices ]; % by default, use all slices
    end

    % if user didn't specify a range of frames to segment, figure it out
    % given the image names
    if ~isfield(p,'segRange')
        D = Dphase1;
        [s,I] = sort({D.name}');
        D = D(I);
        numpos= findstr(D(1).name, '.tif')-3;
        imageNameStrings = char(s);
        p.segRange = str2num(imageNameStrings(:,numpos:numpos+2))';
        clear s I numpos imageNameStrings;
    end
    
%--------------------------------------------------------------------------
% BRIGHTFIELD
elseif (strcmp(p.method,'brightfield')==1)
    Dphase1single   = dir([p.imageDir, '*-p-01-*run-01*.tif']);
    DphaseAllsingle = dir([p.imageDir, '*-p-*run-01*.tif']);
    DphaseAll = dir([p.imageDir, '*-p-*.tif']);
    if ~isfield(p,'numphaseslices')
        if isempty(DphaseAll)
            error(['Can''t find any images in directory ' p.imageDir]);
        else
            p.numphaseslices = int8(ceil(length(DphaseAllsingle)/length(Dphase1single)));
            p.numimagerepetitions = int8(ceil(length(DphaseAll)/length(DphaseAllsingle)));
            disp(['You appear to have ' num2str(p.numphaseslices) ' brightfield images per frame.']);
            disp(['You appear to have ' num2str(p.numimagerepetitions) ' brightfield images per layer per frame.']);
        end
    end

    if ~isfield(p,'slices')
        p.slices = [ 1:p.numphaseslices ]; % by default, use all slices
    end
    if ~any(length(p.slices)==[2 4 8 16 32 64 128])  % if number of slices not a power of 2
        disp('number of slices must be a power of 2! arbitrarily set to 32. pray that it works ...')
        p.slices=[1:32];
    end
    if ~isfield(p,'numimagesperlayer')
        p.numimagesperlayer = [ 1:p.numimagerepetitions ]; % by default, use all slices
    end

    % if user didn't specify a range of frames to segment, figure it out
    % given the image names
    if ~isfield(p,'segRange')
        D = Dphase1single;
        [s,I] = sort({D.name}');
        D = D(I);
        numpos= findstr(D(1).name, '.tif')-3;
        imageNameStrings = char(s);
        p.segRange = str2num(imageNameStrings(:,numpos:numpos+2))';
        clear s I numpos imageNameStrings;
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Display settings that will be used
disp ('using schnitzcells parameter structure:');
disp (p);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Loop over frames : load images and do segmentation
%--------------------------------------------------------------------------
for i= p.segRange
    
    %----------------------------------------------------------------------
    % PHASECONTRAST
    if (strcmp(p.method,'phasecontrast')==1)
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
        %----------------------------------------------------------------------
        % Prepare image and do segmentation
        % ph3   : the phase contrast images of this frame in full size
        % LNsub : the segmented image of this frame in size that contains segmented cells
        % rect  : transformation required to reconstruct full size image from smaller ones

        %if (isfield(p,'medium') & strcmp(p.medium,'rich')==1)        % special
        %                           pre-treatment in rich medium (maybe to come)
        %     imageToSegment = NW_prepareImagesRich(ph3,p.slices);
        %     disp('rich medium')
        %else
        imageToSegment = PN_prepareImages(ph3,p.slices);  
        %end
   
    %----------------------------------------------------------------------
    %BRIGHTFIELD
    elseif (strcmp(p.method,'brightfield')==1)
        % Load images into X and average over repetitions in each layer
        % (per each frame)
        % load averaged layer-images
        if (p.quickMode==1)
            disp('quickMode on. Directly extracting averaged images.')
           
            Dframe = dir([p.imageDir p.NW_brightfieldsaveDir p.movieName '*-p-*-averaged-' ,str3(i), '.tif']); %list of all brightfield layers of this frame (i)
            if isempty(Dframe)
               disp('No averaged pictures found. Will quit...')
               return
            end
        disp('**********************************************')    
        disp(['Up to now no possibility to adjust number of layers/slices used. Will use every layer (', num2str(size(Dframe,1)), ').']);
        disp('**********************************************')
        for runslice=1:size(Dframe,1)
            imagesub(:,:,runslice)=imread([p.imageDir p.NW_brightfieldsaveDir Dframe(runslice).name]);
        end
        
        % load repetitions for each frame and layer and average
        else 
            for runslice=1:p.numphaseslices
                Dframeperlayer = dir([p.imageDir p.movieName '*-p-', num2str(runslice,'%02d'), '*' ,str3(i), '.tif']); %list of all brightfield images of this frame (i) and layer (runslice)

                layername = Dframeperlayer(1).name; % first filename
                fstr = strfind(layername,'-run-');
                for runrep = 1:p.numimagerepetitions
                    layername([fstr+5,fstr+6]) = num2str(runrep,'%02d'); %loads all repetitions of layers (but during averaging less can be used)
                    layerall(:,:,runrep) = imread([p.imageDir,layername]);
                end
                disp(['processing layer (slice) ',num2str(runslice), ' of frame ' num2str(i),'...']);

                % --------------------------------------------------------------------
                % Prepare image: average over repetitions
                layermean = uint16(mean( layerall(:, :, p.numimagesperlayer), 3));
                % save averaged images
                imwrite(uint16(layermean),[p.imageDir, p.NW_brightfieldsaveDir, p.movieName, '-p-', num2str(runslice,'%02d'), '-averaged-',str3(i), '.tif'],'tiff','compression','none'); %blubb change!!!!
                
                % --------------------------------------------------------------------
                % Prepare image: input for correlation analysis
                imageall(:,:,runslice)=layermean;  % matrix with averaged image for each layer
                clear layermean  layerall;
            end
            imagesub=imageall(:,:,p.slices);
            clear imageall ;
        end    
        %----------------------------------------------------------------------
        
        
        % decorrelate images, input: 3dim matrix with image of all layers (imagesub),
        % current frame (i), save directory
        imageToSegment=NW_prepareImages_French_correl_average_frame(p,imagesub,i,[p.imageDir p.NW_correlatedsaveDir]);
        %imageToSegment = PN_prepareImages(ph3,p.slices);  
        
        
        %---------------------------------------------------------------------
        % ph3   : the phase contrast images of this frame in full size
        % LNsub : the segmented image of this frame in size that contains segmented cells
        % rect  : transformation required to reconstruct full size image from smaller ones

    end
    %----------------------------------------------------------------------
    
    %disp(['starting segmentation of frame ',num2str(i), ' Might take forever...']); %blubb
    %return %BLUBB . BREAK FOR TEST PURPOSE ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    saveDirectory = [p.segmentationDir p.PN_saveDir 'seg' str3(i) filesep];
    [status,msg,id] = mkdir(saveDirectory);
    inputsOfSegmentation = {'rangeFiltSize',p.rangeFiltSize,'maskMargin',p.rangeFiltSize,...
        'LoG_Smoothing',p.LoG_Smoothing,'minCellArea',p.minCellArea,...
        'GaussianFilter',p.GaussianFilter,'minDepth',p.minDepth,...
        'neckDepth',p.neckDepth,'saveSteps',p.saveSteps,'saveDir',saveDirectory};
    
    [phsub,LNsub,rect]= NW_segphase_diffStrains(imageToSegment,inputsOfSegmentation{:});                    %%%the real job is done here

    
    %----------------------------------------------------------------------
    % Prepare segmentation and fluorescence data for saving
    if isempty(LNsub) 
        LNsub = [0];          %if no cells segmented, add this so that we have an image
    end

    [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings([p.imageDir,Dframe(1).name],'p'); % extract timestamp of first phase image
    clear TST*

   % phaseFullSize = size( ph3(:,:,1) ); % size of full image
    phaseFullSize = size( imageToSegment(:,:,1) ); % size of full image BLUBB!!!

    savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'''];
    
   % ******************************************************
   % here the fluorescence data was added in former versions
   % ******************************************************
    
    
    % Save segmentation file
    if isfield(p,'PN_saveDir')
        Lname = [p.PN_saveDir, p.movieName, 'seg', str3(i)];
    else
        Lname = [p.movieName, 'seg', str3(i)];
    end
    eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
    disp(['saved file ',p.segmentationDir,Lname]);

    % Save extra seg.png image in \png\ folder
    if isfield(p,'PN_saveDir')
        slices = num2str(p.slices);
        slices(slices==' ') = [];
        png_filename = [png_dir 'seg' str3(i) p.PN_saveDir(7:end-1) '.png'];
        % MW TODO check whether "p.CurrentFrameApprovedFlag = 0;" needs to be
        % added to remove green circles here? 2015/01 
        imwrite(PN_imshowlabel(p,LNsub,[],[],[],'phaseImage',phsub), png_filename, 'png','bitdepth',8);             
    end
end
%-------------------------------------------------------------------------------


