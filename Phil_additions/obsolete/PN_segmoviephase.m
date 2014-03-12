function p = PN_segmoviephase(p,varargin)


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking. Use inputParser in
% recent matlab releases.
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
if ~existfield(p,'cutCellsWidth')                         %minimum neck width to cut a too long cell
    p.cutCellsWidth = 3.5;
end
%saving images
if ~existfield(p,'saveSteps')                             %indicate if you want to save intermediate images
    p.saveSteps = true;
end
if ~existfield(p,'PN_saveDir') & p.saveSteps              %subfolder of p.segmentationDir where image treatment steps are saved
    p.PN_saveDir = ['param' '_Marg' num2str(p.maskMargin) '_LoG' num2str(p.LoG_Smoothing) '_Area' num2str(p.minCellArea) '_Depth' num2str(p.minDepth) '_Width' num2str(p.cutCellsWidth) filesep];
    [message errmsg] = sprintf(['Image saving folder automatically generated: ' p.PN_saveDir]);
    disp(message);
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



%-------------------------------------------------------------------------------
% Figure out which images exist
%--------------------------------------------------------------------------

% figure out the number of phase image layers (slices) per position capture
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
% Display settings that will be used
disp ('using schnitzcells parameter structure:');
disp (p);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Loop over frames : load images and do segmentation
%--------------------------------------------------------------------------
for i= p.segRange
    
    
    %----------------------------------------------------------------------
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
    % Perpare image and do segmentation
    % ph3   : the phase contrast images of this frame in full size
    % LNsub : the segmented image of this frame in size that contains segmented cells
    % rect  : transformation required to reconstruct full size image from smaller ones
    imageToSegment = PN_prepareImages(ph3,p.slices);
    
    saveDirectory = [p.segmentationDir p.PN_saveDir 'seg' str3(i) filesep];
    [status,msg,id] = mkdir(saveDirectory);
    inputsOfSegmentation = {'rangeFiltSize',p.rangeFiltSize,'maskMargin',p.rangeFiltSize,...
        'LoG_Smoothing',p.LoG_Smoothing,'minCellArea',p.minCellArea,...
        'GaussianFilter',p.GaussianFilter,'minDepth',p.minDepth,...
        'cutCellsWidth',p.cutCellsWidth,'saveSteps',p.saveSteps,'saveDir',saveDirectory};
    
    [phsub,LNsub,rect]= PN_segphase(imageToSegment,inputsOfSegmentation{:});                    %%%the real job is done here

    
    %----------------------------------------------------------------------
    % Prepare segmentation and fluorescence data for saving
    if isempty(LNsub) 
        LNsub = [0];          %if no cells segmented, add this so that we have an image
    end

    [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings([p.imageDir,Dframe(1).name],'p'); % extract timestamp of first phase image
    clear TST*

    phaseFullSize = size( ph3(:,:,1) ); % size of full image

    savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'''];

    % Prepare fluor data for saving 
    yname = [p.imageDir,p.movieName,'-y-',str3(i),'.tif'];
    if exist(yname)==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error
        disp('found Fluor image');
        [yreg, yshift, yback, ybinning] = quicknoreg(LNsub,yname,rect,0,phaseFullSize); 
        [exptystr, gainy, expty] = imsettings(yname);
        savelist=[savelist,',''yreg'',''expty'',''gainy'',''yback'',''ybinning'''];
    end
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
        imwrite(DJK_imshowlabel(LNsub,'phaseImage',phsub), png_filename, 'png','bitdepth',8); 
    end
end
%-------------------------------------------------------------------------------


