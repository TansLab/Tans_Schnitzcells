% function DJK_correctFluorImage_anycolor(p, flatfield, shading, replace, varargin)
% CHANGED BY NOREEN
%
% loads fluorimage which is not yet corrected with extra rescale correction
% (for CFP) and neither exists correciton in the .mat data (reg, back,..).
% Does not matter because .mat data will be overwritten!
%
% DJK_correctFluorImage takes original fluor images, corrects them and saves results as:
%  
% example for 'y' fluorescence images:
% * tiff file in '\analysis\fluor_y\y2\' (shading corrected, original size)
% * tiff file in '\analysis\fluor_y\y3\' (shading corrected and shifted, phase contrast size) 
% * tiff file in '\analysis\fluor_y\y4\' (shading corrected and deconvolved, original size)
% * tiff file in '\analysis\fluor_y\y5\' (shading corrected, deconvolved and shifted, phase contrast size) 
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
% 'p'               Note: p.cropLeftTop & p.cropRightBottom must be set
% flatfield         The flatfield matrix (is a misnomer and) contains the
%                   background from microscope settings w. no light to camera 
%                   (historic reasons).
% shading           The shading matrix contains the shading correction image.
% replace           Matrix to indicate dead pixels, 0 means "use this pixel",
%                   1 means "this pixel is dead". Currently _NOT_ in use.
%
% OPTIONAL ARGUMENTS:
% 'fluorcolor'    Fluor color that will be investigated. This can be
%                 'fluor1', 'fluor2', 'fluor3'. The colors associated with
%                 these expressions are stored in 'p'. If no color is
%                 given, 'fluor1'will be chosen.
% 'manualRange'   These frames will be treated
% 'DJK_saveDir'   Folder where tiff files are saved
%                 (default: '\analysis\fluor_[fluorcolor]\')
% 'fluorShift'    Shifts the fluorimage this many pixels
%                 default: [1 -4]
% 'TIFFonly' = 0  default: normal processing: creates TIFF & Fluor using seg files
%            = 1  uses no seg files, and only creates TIFF files. If seg
%                 files do not exist manualRange needs to be given
% 'deconv_func'
% 'rescaleCorrection'  Rescales fluorescence image by an additional factor
%                 (apart from the factor due to different binning) and 
%                 recenters image (according to central image coordinates). In
%                 theory this should not be necessary, but for CFP images
%                 correction is needed (error in optics? etc) . Take around 
%                 XXXX for CFP.
%                 Default: =1
%                 Note that the value of the binning parameter is stored in
%                 the segfile.
% 'minimalMode'       =1: only calculates values, that are actually used
%                     later (Y5xxx,Y6_mean,(+ElowitzStyle) etc).
%                     =0: calculates everything (default: =0)  (NW 2012/04)
%
% IMPLICIT ARGUMENTS
%
% Some arguments used in this function are loaded from the seg file.

function DJK_correctFluorImage_anycolor(p, flatfield, shading, replace, varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_correctFluorImage_anycolor';

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
    % fluorcolor is special case, because the CONTENT of p.fluor1/2/3 must be written to
    % p.fluorcolor and not fluor1/2/3 literally. (NW 11/12/08)
    if strcmp(fieldName,'fluorcolor')==1 
        p.(fieldName)=p.(varargin{i+1});
        disp('--------------------------------------------------------------');
        disp(['Using ' varargin{i+1} ' (' p.(fieldName) ') as fluorescence color.']);
    else
        p.(fieldName) = varargin{i+1};
    end;
  end
end
%--------------------------------------------------------------------------

% If user has set fluorcolor to value 'none', there are no fluor colors to
% correct, and this function will not be executed. MW
if (p.fluorcolor == 'none')
    disp('WARNING, fluorcolor was ''none'', so I WILL NOT CORRECT FLUOR IMAGES.');
    return;
end

%--------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
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

% Choose fluorescence color if not specified
if ~isfield(p,'fluorcolor')
    p.fluorcolor=p.fluor1;
    disp('--------------------------------------------------------------');
    disp(['No fluorescence color specified. Will take fluor1 (' p.fluor1 ').']);
end

% If explicit DJK_saveDir is not given, define it
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'fluor_' p.fluorcolor filesep];
end
% make sure every directory field has a trailing filesep
if (p.DJK_saveDir(length(p.DJK_saveDir)) ~= filesep)
  p.DJK_saveDir = [p.DJK_saveDir filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
fluor2_directory = [p.DJK_saveDir p.fluorcolor '2' filesep];
fluor3_directory = [p.DJK_saveDir p.fluorcolor '3' filesep];
fluor4_directory = [p.DJK_saveDir p.fluorcolor '4' filesep];
fluor5_directory = [p.DJK_saveDir p.fluorcolor '5' filesep];
for i = [1:5]
  switch i,
    case 1, directory = p.DJK_saveDir; 
    case 2, directory = fluor2_directory; 
    case 3, directory = fluor3_directory; 
    case 4, directory = fluor4_directory; 
    case 5, directory = fluor5_directory; 
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
  p.cropLeftTop = [];
  disp(['Warning: p.cropLeftTop was not set. Now []']);
end
if ~existfield(p,'cropRightBottom')
  p.cropRightBottom = []; %[1392,1040];
  disp(['Warning: p.cropRightBottom was not set. Now []']);
end

% If extra rescale Correction for Fluor Image does not exist, set to =1
% (probably all cases except for CFP)
if ~existfield(p,'rescaleCorrection')
    p.rescaleCorrection=1;
end

if ~existfield(p,'minimalMode') % (NW 2012/04)
  p.minimalMode = 0;
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
  disp(['TIFFonly is on. Not using seg files. Will use binning = 2. If error, add manualRange!']);
  binning = 2; % if p.TIFFonly assume that binning is 2
end

disp(['-------------------------------------------------']);
disp(['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
disp(['-------------------------------------------------']);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Prepare flatfield & shading & replace
%--------------------------------------------------------------------------
% rectCrop are coordinates of crop within phaseFullSize
if ~isempty(p.cropLeftTop) & ~isempty(p.cropRightBottom)
    rectCrop = [p.cropLeftTop(2), p.cropLeftTop(1), p.cropRightBottom(2), p.cropRightBottom(1)]; % [top, left, bottom, right]

    % Get correct subset of flatfield & shading
    flatfield_crop  = double( flatfield( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );
    shading_crop    = double(   shading( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );
    replace_crop    = double(   replace( rectCrop(1):rectCrop(3), rectCrop(2):rectCrop(4) ) );
else
    % Do not crop
    rectCrop = [];
    flatfield_crop  = double( flatfield );
    shading_crop    = double(   shading );
    replace_crop    = double(   replace );
end
shading_mean    = mean2(shading);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% check if fluorcolor is non-existent ('none')
%--------------------------------------------------------------------------
if strcmp(p.fluorcolor,'none')==1
    disp('Fluorescence color is non-existent (set to ''none'').');
    optimalShift=[];
    return
else
    %--------------------------------------------------------------------------
    % Generate general variable names (e.g. yreg, ybinning, etc)
    %--------------------------------------------------------------------------
    reg=genvarname([p.fluorcolor 'reg']);
    binning=genvarname([p.fluorcolor 'binning']);
    back=genvarname([p.fluorcolor 'back']);
    gain=genvarname(['gain' p.fluorcolor]);
    expt=genvarname(['expt' p.fluorcolor]);
   
    %--------------------------------------------------------------------------
    % LOOP OVER SEG FILES AND NORMALIZATION OF FLUOR
    %--------------------------------------------------------------------------
    % loop over frames
    for frameNum = p.manualRange
        % load complete fluor image (ALSO THIS WHEN p.TIFFonly)
      fluorname= [p.imageDir,p.movieName,'-' p.fluorcolor '-',str3(frameNum),'.tif'];
      if exist(fluorname)==2
        fluorimage = imread(fluorname);
      else
        disp([' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-' p.fluorcolor '-' str3(frameNum) '.tif in ' p.imageDir]);
        continue;
      end

      % load segmentation file
      if ~p.TIFFonly
        filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];
        % The variables used to be cleared before reassignment but this is not possible with the
        % general names any more. Thus assignment of empty matrices. (NW
        % 11/12/08) <-> think you can -MW
        %eval([reg '=[]; ' binning '=[]; ' back '=[]; ' gain '=[]; ' expt '=[]; ' ]);
        
        % Clear and load parameters used for this frame
        command = ['clear ' reg '; clear ' binning '; clear ' back '; clear ' gain '; clear ' expt '; ']; % MW 2015/04
        eval(command); % MW 2015/04
        clear rect phaseFullSize Lc LNsub;
        load(filename);
            % IMPORTANT NOTE: the segmentation file also stores important 
            % parameters used below, such as e.g. binning.
            % Parameters that are used below are cleared in the above
            % command, where the names of these parameters are stored in 
            % the parameters in the list above.
            % (E.g. the parameter "binning" is a string that holds the name 
            % of the parameter used for actual binning; i.e. this is a 
            % string called ybinning.) MW 2015/04
        disp(['Loaded ' filename]);
 
        disp([' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);
        if ~exist('Lc')      
          disp(['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
          Lc = LNsub;
        end
      end
     
      % still need to resize fluorimage
      eval(['fluorimage = imresize_old(fluorimage,' binning ',''nearest'');']);
      disp(['NOTE: Resizing fluor image by factor ' num2str(eval(binning)) '- MW']);    
       
      %------------------------------------------------------------------------
      % perform manual rescale Correction (maybe unnecessary in this
      % function?) NW 2012-02-27
      %------------------------------------------------------------------------
      if (length(p.rescaleCorrection==1) & p.rescaleCorrection~=1) | length(p.rescaleCorrection)>1
          centerfluor=NW_rescalecenterimage_affine(fluorimage,p.rescaleCorrection);
          fluorimage=centerfluor;
          clear centerfluor;              
      end
      %------------------------------------------------------------------------          
   
      %------------------------------------------------------------------------
      % GET OLD [COLOR]REG & [COLOR]BACK DATA AGAIN
      %------------------------------------------------------------------------
      testreg=eval(['isempty(' reg ')']);
      if (~p.TIFFonly & testreg~=1)
        % Get sizes correct
        phaseCropSize = phaseFullSize; % seg file thinks phaseFullSize is full size, but might be crop size
        phaseFullSize = size(flatfield); % can't do it any other way right now, so use flatfield size as phaseFullSize
        % rect are coordinates of subset within crop
        % rectCrop are coordinates of crop within phaseFullSize
        % rectCropSubset are coordinates of subset within phaseFullSize
        % rectCropSubset = [rectCrop(1)+rect(1), rectCrop(2)+rect(2), rectCrop(1)+rect(3), rectCrop(2)+rect(4)] % [top, left, bottom, right]

        % Background of [color]reg zoals Elowitz doet: median of complete crop fluor image except box around cells
        temp = fluorimage;
        temp(rect(1):rect(3), rect(2):rect(4)) = 0;
        fluorimageVect = temp(temp>0);
        fluorimageVect=double( fluorimageVect);%DE
        eval([back '= uint16( median((fluorimageVect)) );']);

        % Alternative Background: median of fluor within subset (close to cells), but with cells disregarded
        eval([reg '= fluorimage(rect(1):rect(3), rect(2):rect(4) );']);
        eval(['fluorregVect = ' reg '(Lc==0);']);
        fluorregVect=double(fluorregVect);%DE
        fluorbackAlt = uint16( median(fluorregVect) );
      end
      %------------------------------------------------------------------------
      
       %------------------------------------------------------------------------
      % FLATFIELD & SHADING CORRECT
      %------------------------------------------------------------------------
      % Flatfield and Shading correction (ALSO THIS WHEN p.TIFFonly)
      fluor2image = double(fluorimage);
      fluor2image = fluor2image-flatfield_crop; % TODO MW, flatfield should have dependency on illumination time.
      fluor2image = shading_mean.*fluor2image./shading_crop; % when dividing by shading, multiply by mean(shading), so average stays the same
      fluor2image = uint16(fluor2image);
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
      fluor2image_shifted = uint16(zeros(size(fluor2image)));
      fluor2image_shifted = DJK_imageShift(fluor2image_shifted, fluor2image, p.fluorShift);
      %------------------------------------------------------------------------

      %------------------------------------------------------------------------
      % DECONVOLUTION
      %------------------------------------------------------------------------
      eval(['fluor2image_binned = imresize_old(fluor2image,1/' binning ',''nearest'');']);
      fluor2image_deconvolved = p.deconv_func(fluor2image_binned);
      eval(['fluor2image_deconvolved = imresize_old(fluor2image_deconvolved,' binning ',''nearest'');']);
      fluor2image_deconvolved_shifted = uint16(zeros(size(fluor2image_deconvolved)));
      fluor2image_deconvolved_shifted = DJK_imageShift(fluor2image_deconvolved_shifted, fluor2image_deconvolved, p.fluorShift);
      %------------------------------------------------------------------------
      
      %------------------------------------------------------------------------
      % GET FLATFIELD & SHADING CORRECTED [COLOR]REG2 & [COLOR]BACK2 DATA
      %------------------------------------------------------------------------
      testreg=eval(['isempty(' reg ')']);
      if (~p.TIFFonly & testreg~=1)
          
        if p.minimalMode~=1
            % Background of reg2 zoals Elowitz doet: median of complete crop fluor image except box around cells
            temp = fluor2image;
            temp(rect(1):rect(3), rect(2):rect(4)) = 0;
            fluorimageVect = temp(temp>0);
            fluorimageVect=double(fluorimageVect);%DE 
            fluorback2 = uint16( median(fluorimageVect) );

            % Alternative Background: median of fluorcolor within subset (close to cells), but with cells disregarded
            fluorreg2 = fluor2image(rect(1):rect(3), rect(2):rect(4) );
            fluorreg2Vect = fluorreg2(Lc==0);
            fluorreg2Vect=double(fluorreg2Vect);
            fluorback2Alt = uint16( median(fluorreg2Vect) ); 
        end

        % Background of fluorreg3 zoals Elowitz doet: median of complete crop fluor image except box around cells
        temp = fluor2image_shifted;
        temp(rect(1):rect(3), rect(2):rect(4)) = 0;
        fluorimageVect = temp(temp>0);
        fluorimageVect=double(fluorimageVect);%DE
        fluorback3 = uint16( median(fluorimageVect) );

        % Get correct subset
        fluorreg3 = fluor2image_shifted(rect(1):rect(3), rect(2):rect(4) );
        fluorreg4 = fluor2image_deconvolved(rect(1):rect(3), rect(2):rect(4) );
        fluorreg5 = fluor2image_deconvolved_shifted(rect(1):rect(3), rect(2):rect(4) );
      end
      %------------------------------------------------------------------------
      
      
      %----------------------------------------------------------------------
      % CREATE CORRECT (COLOR SPECIFIC) VARIABLE NAMES AND ASSOCIATE DATA (maybe do at
      % beginning) NW 11/12/08
      %----------------------------------------------------------------------
      shift=genvarname([p.fluorcolor 'shift']);
      reg5=genvarname([p.fluorcolor 'reg5']); eval([reg5 '=fluorreg5;']);
      back3=genvarname([p.fluorcolor 'back3']); eval([back3 '=fluorback3;']);    
      if p.minimalMode~=1
          reg2=genvarname([p.fluorcolor 'reg2']); eval([reg2 '=fluorreg2;']);
          reg3=genvarname([p.fluorcolor 'reg3']); eval([reg3 '=fluorreg3;']);
          reg4=genvarname([p.fluorcolor 'reg4']); eval([reg4 '=fluorreg4;']);
          back2=genvarname([p.fluorcolor 'back2']); eval([back2 '=fluorback2;']);
         backAlt=genvarname([p.fluorcolor 'backAlt']); eval([backAlt '=fluorbackAlt;']);
          back2Alt=genvarname([p.fluorcolor 'back2Alt']); eval([back2Alt '=fluorback2Alt;']);
      end
      f2image_binned=genvarname([p.fluorcolor '2image_binned']); eval([f2image_binned '=fluor2image_binned;']);
      f2image_shifted=genvarname([p.fluorcolor '2image_shifted']); eval([f2image_shifted '=fluor2image_shifted;']);
      f2image_deconvolved_binned=genvarname([p.fluorcolor '2image_deconvolved_binned']);% eval([f2image_deconvolved_binned '=fluor2image_deconvolved_binned;']);
      f2image_deconvolved_shifted=genvarname([p.fluorcolor '2image_deconvolved_shifted']); eval([f2image_deconvolved_shifted '=fluor2image_deconvolved_shifted;']);
      
   
      %----------------------------------------------------------------------
      % SAVE TO MAT FILE
      %----------------------------------------------------------------------
 
      testreg=eval(['isempty(' reg ')']);
      if (~p.TIFFonly & testreg~=1)
        eval([shift '= p.fluorShift;']);
        if p.minimalMode~=1
            savelist = ['''rectCrop'',''rect'',''phaseFullSize'',''phaseCropSize'',''' gain ''',''' binning ''',''' expt ''',''' shift ''',''' reg ''',''' reg2 ''',''' reg3 ''',''' reg4 ''',''' reg5 ''',''' back ''',''' back2 ''',''' back3 ''',''' backAlt ''',''' back2Alt ''''];
        else
             savelist = ['''rectCrop'',''rect'',''phaseFullSize'',''phaseCropSize'',''' gain ''',''' binning ''',''' expt ''',''' shift ''',''' reg ''',''' reg5 ''',''' back ''',''' back3  ''''];
        end
        filename = [p.tracksDir, p.movieName, 'Fluor_', p.fluorcolor '_', str3(frameNum)]; 
        eval(['save(''' filename ''',' savelist ');']);
        disp(['       -> saved ' [p.movieName, 'Fluor_', p.fluorcolor '_', str3(frameNum)] ' in ' p.tracksDir]);
      end;
      %----------------------------------------------------------------------
      
      
      %----------------------------------------------------------------------
      % SAVE TIFF IMAGES OF FLUOR2
      %----------------------------------------------------------------------
        eval([f2image_binned '= imresize_old(fluor2image ,1/' binning ',''nearest'');']) %repeated?
        str=['imwrite(' f2image_binned ', ''' [ fluor2_directory  p.movieName '-' p.fluorcolor '2-' str3(frameNum) '.tif' ] ''' , ''TIFF'', ''Compression'', ''none'');'];
        eval(str);
        disp(['       -> saved ' [p.movieName, '-', p.fluorcolor, '2-', str3(frameNum) '.tif'] ' in ' fluor2_directory]);
        str=['imwrite(' f2image_shifted ', ''' [ fluor3_directory  p.movieName '-' p.fluorcolor '3-' str3(frameNum) '.tif' ] ''' , ''TIFF'', ''Compression'', ''none'');'];
        eval(str);
        disp(['       -> saved ' [p.movieName, '-', p.fluorcolor, '3-', str3(frameNum) '.tif'] ' in ' fluor3_directory]);
        eval([f2image_deconvolved_binned ' = imresize_old( fluor2image_deconvolved ,1/' binning ',''nearest'');']);
        str=['imwrite(' f2image_deconvolved_binned ', ''' [ fluor4_directory  p.movieName '-' p.fluorcolor '4-' str3(frameNum) '.tif' ] ''' , ''TIFF'', ''Compression'', ''none'');'];
        eval(str);
        disp(['       -> saved ' [p.movieName, '-', p.fluorcolor, '4-', str3(frameNum) '.tif'] ' in ' fluor4_directory]);
        str=['imwrite(' f2image_deconvolved_shifted ', ''' [ fluor5_directory  p.movieName '-' p.fluorcolor '5-' str3(frameNum) '.tif' ] ''' , ''TIFF'', ''Compression'', ''none'');'];
        eval(str);
        disp(['       -> saved ' [p.movieName, '-', p.fluorcolor, '5-', str3(frameNum) '.tif'] ' in ' fluor5_directory]);
       
      
        %----------------------------------------------------------------------
    end
    %--------------------------------------------------------------------------

    disp(['-------------------------------------------------']);

  

end

