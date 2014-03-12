function DJK_cropImages(p, cropRange, leftTop, rightBottom, varargin) 
% DJK_cropImages creates crops of all images. Coordinates given must be even (in
% order for the fluor images to be cropped correctly)
% 
%   'cropRange'       Defines a specific range of frame numbers to extract. 
%                     By default, new movie will contain as many frames as 
%                     the given original movie.
%
%   'leftTop'         Coordinate of left top of crop (must be even).
%
%   'rightBottom'     Coordinate of right bottom of crop (must be even).
%
%   'cropName'        New name for cropped images (standard: posXcrop)
%


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4;
functionName = 'DJK_cropImages';

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

% check that coordinates are even
%--------------------------------------------------------------------------
if (~mod(leftTop(1),2) || ~mod(leftTop(2),2))
      errorMessage = sprintf('%s\n',['leftTop should contain uneven values: ' num2str(leftTop(1)) ',' num2str(leftTop(2))]);
      error(errorMessage);
end
if (mod(rightBottom(1),2) || mod(rightBottom(2),2)) 
      errorMessage = sprintf('%s\n',['rightBottom should contain even values: ' num2str(rightBottom(1)) ',' num2str(rightBottom(2))]);
      error(errorMessage);
end

%--------------------------------------------------------------------------
% Parse the input arguments
% Override any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(p,'cropName')
  p.cropName = [p.movieName, 'crop'];
end

if ~existfield(p,'cropDir')
  p.cropDir = [p.dateDir, p.cropName, '\'];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Figure out what images are there
%--------------------------------------------------------------------------
% All phase images
DphaseAll = dir([p.imageDir, [p.movieName,'-p-*.tif'] ]);
% All phase -p-2- images
Dphase2   = dir([p.imageDir, [p.movieName,'-p-2-*.tif'] ]);
% All fluor y images
DfluorY = dir([p.imageDir, [p.movieName,'-y-*.tif'] ]);
% All fluor r images
DfluorR = dir([p.imageDir, [p.movieName,'-r-*.tif'] ]);

if isempty(Dphase2)
    disp('You seem to have 1 phase image per frame');
    if ~isempty(DfluorY)
        disp(['You seem to a fluor image for every ' num2str(length(DphaseAll)/length(DfluorY)) ' frames']);
    end 
    if ~isempty(DfluorR)
        disp(['You seem to a fluor2 image for every ' num2str(length(DphaseAll)/length(DfluorY)) ' frames']);
    end 
else 
    disp(['You seem to have ' num2str(length(DphaseAll)/length(Dphase2)) ' phase images per frame']);
    if ~isempty(DfluorY)
        disp(['You seem to a fluor image for every ' num2str(length(Dphase2)/length(DfluorY)) ' frames ']);
    end 
    if ~isempty(DfluorY)
        disp(['You seem to a fluor2 image for every ' num2str(length(Dphase2)/length(DfluorY)) ' frames ']);
    end 
end
if isempty(DfluorY)
    disp('You seem to have no fluor images');
end
if isempty(DfluorR)
    disp('You seem to have no fluor2 images');
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Create Crop folder where images will be saved
%--------------------------------------------------------------------------
% make sure directory field has a trailing filesep
if (p.cropDir( length(p.cropDir) ) ~= filesep)
  p.cropDir = [p.cropDir filesep];
end

% if directories don't exist, create them (will issue warnings if they exist)
if exist(p.cropDir)~=7
  [status,msg,id] = DJK_mkdir(p.cropDir);
  if status == 0
    disp(['Warning: unable to mkdir ' p.cropDir ' : ' msg]);
  end
end

cropDirImageDir = [p.cropDir 'images\'];
if exist(cropDirImageDir)~=7
  [status,msg,id] = DJK_mkdir(cropDirImageDir);
  if status == 0
    disp(['Warning: unable to mkdir ' cropDirImageDir ' : ' msg]);
  end
else 
  disp(['cropDir already exists -> could be overwriting files in: ' cropDirImageDir]);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Perform Cropping
%--------------------------------------------------------------------------
for fr = cropRange, % go over each frame
    DphaseRange = dir([p.imageDir, [p.movieName, '-p-*', str3(fr) ,'.tif'] ]);
    DfluorYRange = dir([p.imageDir, [p.movieName, '-y-*', str3(fr) ,'.tif'] ]);
    DfluorRRange = dir([p.imageDir, [p.movieName, '-r-*', str3(fr) ,'.tif'] ]);
    
    for i = [1:length(DphaseRange)], % go over each phase image of this frame
        % read image
        im_original = imread([p.imageDir DphaseRange(i).name]); 
        % get image info
        im_info = imfinfo([p.imageDir DphaseRange(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        % crop the image
        im_crop = im_original( leftTop(2):rightBottom(2), leftTop(1):rightBottom(1));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(DphaseRange(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written: ' im_crop_filename]);
        %disp([ num2str(leftTop(2)) ':' num2str(rightBottom(2)) '-' num2str(leftTop(1)) ':' num2str(rightBottom(1))]);
    end
    
    for i = [1:length(DfluorYRange)],  % go over each fluor image of this frame
        % read image
        im_original = imread([p.imageDir DfluorYRange(i).name]); 
        % get image info
        im_info = imfinfo([p.imageDir DfluorYRange(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(DfluorYRange(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written: ' im_crop_filename]);
        %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
    for i = [1:length(DfluorRRange)],  % go over each fluor image of this frame
        % read image
        im_original = imread([p.imageDir DfluorRRange(i).name]); 
        % get image info
        im_info = imfinfo([p.imageDir DfluorRRange(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(DfluorRRange(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written: ' im_crop_filename]);
        %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
end
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [status,msg,id] = DJK_mkdir(absoluteDir);
% This is a wrapper around mkdir, because mkdir with one argument 
% behaves differently in matlab 6.5 vs. matlab 7.0.  This behaves
% like matlab 7.0 mkdir when called with one argument, and it assumes 
% the one argument is an absolute directory.  Works on Win & Unix.
v = version;
if (v(1)=='7')
  [status, msg, id] = mkdir(absoluteDir);
  return
end
% below, we assume we're in version matlab 6.*
if isunix
  [status, msg, id] = mkdir('/',absoluteDir);
  return
end
if ispc
  [status, msg, id] = mkdir(absoluteDir(1:2),absoluteDir(3:end));
  return
end
error('mymkdir cant figure out what OS you have!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
