function PN_copySegFiles(p,varargin)
% function PN_copySegFiles(p,varargin)
%
% (comments by MW 2015/06)
% Copies segmentation files from subdirectory with annotated parameter
% values to main dir with segmentation files.
%
% overwrite     defualt=1, will overwrite old segfiles. if overwrite=0, old
%               segfiles won't be overwritten.

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

%%%%%%%default values. Have to be similar to the ones in PN_segmoviephase.
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
if ~existfield(p,'neckDepth')                         %minimum neck width to cut a too long cell
    p.neckDepth = 2;
end
if ~isfield(p,'segRange')
    error('Field segRange empty : please specify the images number to copy.')
end
if ~existfield(p,'overwrite')                         %minimum neck width to cut a too long cell
    p.overwrite = 0;
end


%--------------------------------------------------------------------------
%copy files to be used in following steps in the segmentation folder
%--------------------------------------------------------------------------

directory = [p.segmentationDir 'param' '_Marg' num2str(p.maskMargin) '_LoG' num2str(p.LoG_Smoothing) '_Area' num2str(p.minCellArea) '_Depth' num2str(p.minDepth) '_Neck' num2str(p.neckDepth) filesep];
    if ~isdir(directory)
        errormessage = sprintf('%s\n%s\n','The following directory does not exist :',directory);
        error(errormessage);
    end

for i= p.segRange   
    filename = [p.movieName 'seg' str3(i) '.mat'];
    if exist([directory filename],'file')==0
        warningmessage = sprintf('%s\n%s\n','The following file does not exist :',[directory filename]);
        warning(warningmessage);
    end
       
    % MW addition to keep already existing segfiles
    if p.overwrite==0 && exist([p.segmentationDir filename],'file')~=0
        disp(['Didn''t copy ' filename ' because old segfile already detected (p.overwrite==0).']);
        continue;
    elseif exist([p.segmentationDir filename],'file')~=0
        disp('WARNING: old segfile is being overwritten.');
    end    
    
    disp(['Copying ' filename '..']);
    try
        copyfile([directory filename],[p.segmentationDir filename]);                % copy *.mat files
    catch 
        disp(['ERROR: File not found: ' filename '; Continuing anyway.']);
    end
    
end



