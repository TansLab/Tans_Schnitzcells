% DJK_correctFluorImage takes original fluor images, corrects them and saves results as:
% * tiff file in '\analysis\fluor\r2\' (shading corrected, original size)
% * tiff file in '\analysis\fluor\r3\' (shading corrected and shifted, phase contrast size) 
% * tiff file in '\analysis\fluor\r4\' (shading corrected and deconvolved, original size)
% * tiff file in '\analysis\fluor\r5\' (shading corrected, deconvolved and shifted, phase contrast size) 
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
% 'DJK_saveDir2'   Folder where tiff files are saved
%                 (default: '\analysis\fluor\')
% 'fluorShift'    Shifts the y image this many pixels
%                 default: [1 -4]
% 'TIFFonly' = 0  default: normal processing: creates TIFF & Fluor using seg files
%            = 1  uses no seg files, and only creates TIFF files. If seg
%                 files do not exist manualRange needs to be given
% 'deconv_func'

function DJK_correctFluorImage_red(p, flatfield, shading, replace, varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_correctFluorImage_red';

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

% If explicit DJK_saveDir2 is not given, define it
if ~existfield(p,'DJK_saveDir2')
  p.DJK_saveDir2 = [p.analysisDir 'fluor2\'];
end
% make sure every directory field has a trailing filesep
if (p.DJK_saveDir2(length(p.DJK_saveDir2)) ~= filesep)
  p.DJK_saveDir2 = [p.DJK_saveDir2 filesep];
end
r2_directory = [p.DJK_saveDir2 'r2\'];
r3_directory = [p.DJK_saveDir2 'r3\'];
r4_directory = [p.DJK_saveDir2 'r4\'];
r5_directory = [p.DJK_saveDir2 'r5\'];
for i = [1:5]
  switch i,
    case 1, directory = p.DJK_saveDir2; 
    case 2, directory = r2_directory; 
    case 3, directory = r3_directory; 
    case 4, directory = r4_directory; 
    case 5, directory = r5_directory; 
  end
  if exist(directory)~=7
    [status,msg,id] = mymkdir([directory]);
    if status == 0
      disp(['Warning: unable to mkdir ' directory ' : ' msg]); return;
    end
  end
end

% If explicit shift is not given, getShift is [0 0]
if ~existfield(p,'fluor2Shift')
  p.fluor2Shift = [0 0];
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
disp([' * fluor2Shift [' num2str(p.fluor2Shift(1)) ',' num2str(p.fluor2Shift(2)) ']']);

if p.TIFFonly
  disp(['-------------------------------------------------']);
  disp(['TIFFonly is on. Not using seg files. Will use rbinning = 2. If error, add manualRange!']);
  rbinning = 2; % if p.TIFFonly assume that binning is 2
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
  rname= [p.imageDir,p.movieName,'-r-',str3(frameNum),'.tif'];
  if exist(rname)==2
    rimage = imread(rname);
  else
    disp([' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-r-' str3(frameNum) '.tif in ' p.imageDir]);
    continue;
  end

  % load segmentation file
  if ~p.TIFFonly
    filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];
    clear rreg rect rbinning rback phaseFullSize gainr exptr Lc LNsub;
    rreg = [];
    load(filename);
    disp([' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);
    if ~exist('Lc')      
      disp(['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
      Lc = LNsub;
    end
  end

  % still need to resize rimage
  rimage = imresize_old(rimage,rbinning,'nearest');
  
  %------------------------------------------------------------------------
  % GET OLD YREG & YBACK DATA AGAIN
  %------------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(rreg)
    % Get sizes correct
    phaseCropSize = phaseFullSize; % seg file thinks phaseFullSize is full size, but might be crop size
    phaseFullSize = size(flatfield); % can't do it any other way right now, so use flatfield size as phaseFullSize
    % rect are coordinates of subset within crop
    % rectCrop are coordinates of crop within phaseFullSize
    % rectCropSubset are coordinates of subset within phaseFullSize
    % rectCropSubset = [rectCrop(1)+rect(1), rectCrop(2)+rect(2), rectCrop(1)+rect(3), rectCrop(2)+rect(4)] % [top, left, bottom, right]

    % Background of yreg zoals Elowitz doet: median of complete crop y image except box around cells
    temp = rimage;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    rimageVect = temp(temp>0);
    rback = uint16( median(rimageVect) );

    % Alternative Background: median of y within subset (close to cells), but with cells disregarded
    rreg = rimage(rect(1):rect(3), rect(2):rect(4) );
    rregVect = rreg(Lc==0);
    rbackAlt = uint16( median(rregVect) );
  end
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % FLATFIELD & SHADING CORRECT
  %------------------------------------------------------------------------
  % Flatfield and Shading correction (ALSO THIS WHEN p.TIFFonly)
  r2image = double(rimage);
  r2image = r2image-flatfield_crop;
  r2image = shading_mean.*r2image./shading_crop; % when dividing by shading, multiply by mean(shading), so average stays the same
  r2image = uint16(r2image);
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
  r2image_shifted = uint16(zeros(size(r2image)));
  r2image_shifted = DJK_imageShift(r2image_shifted, r2image, p.fluor2Shift);
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % DECONVOLUTION
  %------------------------------------------------------------------------
  r2image_binned = imresize_old(r2image,1/rbinning,'nearest');
  r2image_deconvolved = p.deconv_func(r2image_binned);
  r2image_deconvolved = imresize_old(r2image_deconvolved,rbinning,'nearest');
  r2image_deconvolved_shifted = uint16(zeros(size(r2image_deconvolved)));
  r2image_deconvolved_shifted = DJK_imageShift(r2image_deconvolved_shifted, r2image_deconvolved, p.fluor2Shift);
  %------------------------------------------------------------------------
  
  %------------------------------------------------------------------------
  % GET FLATFIELD & SHADING CORRECTED RREG2 & RBACK2 DATA
  %------------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(rreg)

    % Background of rreg2 zoals Elowitz doet: median of complete crop r image except box around cells
    temp = r2image;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    rimageVect = temp(temp>0);
    rback2 = uint16( median(rimageVect) );

    % Alternative Background: median of r within subset (close to cells), but with cells disregarded
    rreg2 = r2image(rect(1):rect(3), rect(2):rect(4) );
    rreg2Vect = rreg2(Lc==0);
    rback2Alt = uint16( median(rreg2Vect) ); 

    % Background of yreg3 zoals Elowitz doet: median of complete crop r image except box around cells
    temp = r2image_shifted;
    temp(rect(1):rect(3), rect(2):rect(4)) = 0;
    rimageVect = temp(temp>0);
    rback3 = uint16( median(rimageVect) );

    % Get correct subset
    rreg3 = r2image_shifted(rect(1):rect(3), rect(2):rect(4) );
    rreg4 = r2image_deconvolved(rect(1):rect(3), rect(2):rect(4) );
    rreg5 = r2image_deconvolved_shifted(rect(1):rect(3), rect(2):rect(4) );
  end
  %------------------------------------------------------------------------

  %----------------------------------------------------------------------
  % SAVE TO MAT FILE
  %----------------------------------------------------------------------
  if ~p.TIFFonly & ~isempty(rreg)
    rshift = p.fluor2Shift;
    savelist = ['''rectCrop'',''rect'',''phaseFullSize'',''phaseCropSize'',''gainr'',''rbinning'',''exptr'',''rshift'',''rreg'',''rreg2'',''rreg3'',''rreg4'',''rreg5'',''rback'',''rback2'',''rback3'',''rbackAlt'',''rback2Alt'''];
    filename = [p.tracksDir, p.movieName, 'Fluor2_', str3(frameNum)];
    eval(['save(''' filename ''',' savelist ');']);
    disp(['       -> saved ' [p.movieName, 'Fluor2_', str3(frameNum)] ' in ' p.tracksDir]);
  end
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % SAVE TIFF IMAGES OF R2
  %----------------------------------------------------------------------
    r2image_binned = imresize_old(r2image,1/rbinning,'nearest');
    imwrite(r2image_binned, [r2_directory, p.movieName, '-r2-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-r2-', str3(frameNum) '.tif'] ' in ' r2_directory]);
    imwrite(r2image_shifted, [r3_directory, p.movieName, '-r3-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-r3-', str3(frameNum) '.tif'] ' in ' r3_directory]);
    r2image_deconvolved_binned = imresize_old(r2image_deconvolved,1/rbinning,'nearest');
    imwrite(r2image_deconvolved_binned, [r4_directory, p.movieName, '-r4-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-r4-', str3(frameNum) '.tif'] ' in ' r4_directory]);
    imwrite(r2image_deconvolved_shifted, [r5_directory, p.movieName, '-r5-', str3(frameNum) '.tif'], 'TIFF', 'Compression', 'none');
    disp(['       -> saved ' [p.movieName, '-r5-', str3(frameNum) '.tif'] ' in ' r5_directory]);
  %----------------------------------------------------------------------
end 
%--------------------------------------------------------------------------

disp(['-------------------------------------------------']);
