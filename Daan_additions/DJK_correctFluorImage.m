% DJK_correctFluorImage takes original fluor images, corrects them and saves results as:
% * tiff file in '\analysis\fluor\y2\' (shading corrected, original size)
% * tiff file in '\analysis\fluor\y3\' (shading corrected and shifted, phase contrast size) 
% * tiff file in '\analysis\fluor\y4\' (shading corrected and deconvolved, original size)
% * tiff file in '\analysis\fluor\y5\' (shading corrected, deconvolved and shifted, phase contrast size) 
% * xxxFluorxxx.mat in folder '\data\'
%
% The following corrections are applied: 
% * using flatfield, shading & replace corrects for shading / flatfield /
%   dead pixels
% * when deconv_func argument is provided, the accompanying funtion & PSF is used
%   for deconvolution
% * using fluorShift (determined by DJK_getFluorShift) the fluor image is also shifted
%
% OUTPUT
%
% REQUIRED ARGUMENTS:
% 'p'             Note: p.cropLeftTop & p.cropRightBottom must be set
%
% OPTIONAL ARGUMENTS:
% 'manualRange'   These frames will be treated
% 'DJK_saveDir'   Folder where tiff files are saved
%                 (default: '\analysis\fluor\')
% 'fluorShift'    Shifts the y image this many pixels
%                 default: [1 -4]
% 'TIFFonly' = 0  default: normal processing: creates TIFF & Fluor using seg files
%            = 1  uses no seg files, and only creates TIFF files. If seg
%                 files do not exist manualRange needs to be given
% 'deconv_func'

function DJK_correctFluorImage(p, flatfield, shading, replace, varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_correctFluorImage';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
  error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If explicit manualRange is not given, take all segmentation files
if ~existfield(p,'manualRange')
  % Get directory of existing segmentation files 
  outprefix = [p.movieName 'seg'];
  D = dir([p.segmentationDir, outprefix, '*.mat']);
  [S,I] = sort({D.name}');
  D = D(I);
  numpos= findstr(D(1).name, '.mat')-3;
  
  segNameStrings = char(S);
  p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
end

% If explicit DJK_saveDir is not given, define it
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'fluor\'];
end
% make sure every directory field has a trailing filesep
if (p.DJK_saveDir(length(p.DJK_saveDir)) ~= filesep)
  p.DJK_saveDir = [p.DJK_saveDir filesep];
end
y2_directory = [p.DJK_saveDir 'y2\'];
y3_directory = [p.DJK_saveDir 'y3\'];
y4_directory = [p.DJK_saveDir 'y4\'];
y5_directory = [p.DJK_saveDir 'y5\'];
for i = [1:5]
  switch i,
    case 1, directory = p.DJK_saveDir; 
    case 2, directory = y2_directory; 
    case 3, directory = y3_directory; 
    case 4, directory = y4_directory; 
    case 5, directory = y5_directory; 
  end
  if exist(directory)~=7
    [status,msg,id] = mymkdir([directory]);
    if status == 0
      disp(['Warning: unable to mkdir ' directory ' : ' msg]); return;
    end
  end
end

% If explicit shift is not given, getShift is [0 0]
if ~existfield(p,'fluorShift')
  p.fluorShift = [0 0];
end

% If explicit TIFFonly is not given, TIFFonly is false (0)
if ~existfield(p,'TIFFonly')
  p.TIFFonly = 0;
end

% If explicit deconv_func is not given, deconv_func returns the same images
if ~existfield(p,'deconv_func')
  p.deconv_func = @(im) im;
end

% Check whether cropLeftTop & cropRightBottom are set
if ~existfield(p,'cropLeftTop')
  cropLeftTop = [1,1];
  disp(['Warning: p.cropLeftTop was not set. Now [1,1]']);
end
if ~existfield(p,'cropRightBottom')
  cropRightBottom = [1392,1040];
  disp(['Warning: p.cropRightBottom was not set. Now [1392,1040]']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Let know what is happening
%--------------------------------------------------------------------------
disp(['-------------------------------------------------']);
disp(['Correcting fluor image with:']);
disp([' * flatfield']);
disp([' * shading']);
disp([' * replace']);
disp([' * fluorShift [' num2str(p.fluorShift(1)) ',' num2str(p.fluorShift(2)) ']']);

if p.TIFFonly
  disp(['-------------------------------------------------']);
  disp(['TIFFonly is on. Not using seg files. Will use ybinning = 2. If error, add manualRange!']);
  ybinning = 2; % if p.TIFFonly assume that binning is 2
end

disp(['-------------------------------------------------']);
disp(['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
disp(['-------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Prepare flatfield & shading & replace
%--------------------------------------------------------------------------
% rectCrop are coordinates of crop within phaseFullSize
rectCrop = [p.cropLeftTop(2), p.cropLeftTop(1), p.cropRightBottom(2), p.cropRightBottom(1)]; % [top, left, bottom, right]

% Get correct subset of flatfield & shading
flatfield_crop  = double( flatfield( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );
shading_crop    = double(   shading( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );
replace_crop    = double(   replace( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );

shading_mean    = mean(mean(shading));
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SEG FILES AND NORMALIZATION OF FLUOR
%--------------------------------------------------------------------------
% loop over frames
for frameNum = p.manualRange
  
  % load complete fluor image (ALSO THIS WHEN p.TIFFonly)
  yname= [p.imageDir,p.movieName,'-y-',str3(frameNum),'.tif'];
  if exist(yname)==2
    yimage = imread(yname);
  else
    disp([' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-y-' str3(frameNum) '.tif in ' p.imageDir]);
    continue;
  end

  % load segmentation file
  if ~p.TIFFonly
    filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];
    clear yreg rect ybinning yback phaseFullSize gainy expty Lc LNsub;
    yreg = [];
    load(filename);
    disp([' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);
    if ~exist('Lc')      
      disp(['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
      Lc = LNsub;
    end
  end

  % still need to resize yimage
  yimage = imresize_old(yimage,ybinning,'nearest');
  
  %------------------------------------------------------------------------
  % GET OLD YREG & YBACK DATA AGAIN
  %------------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(yreg)
    % Get sizes correct
    phaseCropSize = phaseFullSize; % seg file thinks phaseFullSize is full size, but might be crop size
    phaseFullSize = size(flatfield); % can't do it any other way right now, so use flatfield size as phaseFullSize
    % rect are coordinates of subset within crop
    % rectCrop are coordinates of crop within phaseFullSize
    % rectCropSubset are coordinates of subset within phaseFullSize
    % rectCropSubset = [rectCrop(1)+rect(1), rectCrop(2)+rect(2), rectCrop(1)+rect(3), rectCrop(2)+rect(4)] % [top, left, bottom, right]

    % Background of yreg zoals Elowitz doet: median of complete crop y image except box around cells
    temp = yimage;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    yimageVect = temp(temp>0);
    yback = uint16( median(yimageVect) );

    % Alternative Background: median of y within subset (close to cells), but with cells disregarded
    yreg = yimage(rect(1):rect(3), rect(2):rect(4) );
    yregVect = yreg(Lc==0);
    ybackAlt = uint16( median(yregVect) );
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FLATFIELD & SHADING CORRECT
  %------------------------------------------------------------------------
  % Flatfield and Shading correction (ALSO THIS WHEN p.TIFFonly)
  y2image = double(yimage);
  y2image = y2image-flatfield_crop;
  y2image = shading_mean.*y2image./shading_crop; % when dividing by shading, multiply by mean(shading), so average stays the same
  y2image = uint16(y2image);
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % REPLACE DEAD PIXELS CORRECT
  %------------------------------------------------------------------------
  % 2DO
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % SHIFT Y IMAGE
  %------------------------------------------------------------------------
  % Shift 
  y2image_shifted = uint16(zeros(size(y2image)));
  y2image_shifted = DJK_imageShift(y2image_shifted, y2image, p.fluorShift);
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % DECONVOLUTION
  %------------------------------------------------------------------------
  y2image_binned = imresize_old(y2image,1/ybinning,'nearest');
  y2image_deconvolved = p.deconv_func(y2image_binned);
  y2image_deconvolved = imresize_old(y2image_deconvolved,ybinning,'nearest');
  y2image_deconvolved_shifted = uint16(zeros(size(y2image_deconvolved)));
  y2image_deconvolved_shifted = DJK_imageShift(y2image_deconvolved_shifted, y2image_deconvolved, p.fluorShift);
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % GET FLATFIELD & SHADING CORRECTED YREG2 & YBACK2 DATA
  %------------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(yreg)

    % Background of yreg2 zoals Elowitz doet: median of complete crop y image except box around cells
    temp = y2image;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    yimageVect = temp(temp>0);
    yback2 = uint16( median(yimageVect) );

    % Alternative Background: median of y within subset (close to cells), but with cells disregarded
    yreg2 = y2image(rect(1):rect(3), rect(2):rect(4) );
    yreg2Vect = yreg2(Lc==0);
    yback2Alt = uint16( median(yreg2Vect) ); 

    % Background of yreg3 zoals Elowitz doet: median of complete crop y image except box around cells
    temp = y2image_shifted;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    yimageVect = temp(temp>0);
    yback3 = uint16( median(yimageVect) );

    % Get correct subset
    yreg3 = y2image_shifted(rect(1):rect(3), rect(2):rect(4) );
    yreg4 = y2image_deconvolved(rect(1):rect(3), rect(2):rect(4) );
    yreg5 = y2image_deconvolved_shifted(rect(1):rect(3), rect(2):rect(4) );
  end
  %------------------------------------------------------------------------

  %----------------------------------------------------------------------
  % SAVE TO MAT FILE
  %----------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(yreg)
    yshift = p.fluorShift;
    savelist = ['''rectCrop'',''rect'',''phaseFullSize'',''phaseCropSize'',''gainy'',''ybinning'',''expty'',''yshift'',''yreg'',''yreg2'',''yreg3'',''yreg4'',''yreg5'',''yback'',''yback2'',''yback3'',''ybackAlt'',''yback2Alt'''];
    filename = [p.tracksDir, p.movieName, 'Fluor', str3(frameNum)];
    eval(['save(''' filename ''',' savelist ');']);
    disp(['       -> saved ' [p.movieName, 'Fluor', str3(frameNum)] ' in ' p.tracksDir]);
  end
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % SAVE TIFF IMAGES OF Y2
  %----------------------------------------------------------------------
    y2image_binned = imresize_old(y2image,1/ybinning,'nearest');
    imwrite(y2image_binned, [y2_directory, p.movieName, '-y2-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-y2-', str3(frameNum) '.tif'] ' in ' y2_directory]);
    imwrite(y2image_shifted, [y3_directory, p.movieName, '-y3-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-y3-', str3(frameNum) '.tif'] ' in ' y3_directory]);
    y2image_deconvolved_binned = imresize_old(y2image_deconvolved,1/ybinning,'nearest');
    imwrite(y2image_deconvolved_binned, [y4_directory, p.movieName, '-y4-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-y4-', str3(frameNum) '.tif'] ' in ' y4_directory]);
    imwrite(y2image_deconvolved_shifted, [y5_directory, p.movieName, '-y5-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-y5-', str3(frameNum) '.tif'] ' in ' y5_directory]);
  %----------------------------------------------------------------------
end 
%--------------------------------------------------------------------------

disp(['-------------------------------------------------']);
