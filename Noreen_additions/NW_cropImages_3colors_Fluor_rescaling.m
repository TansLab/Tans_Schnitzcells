function NW_cropImages_3colors_Fluor_rescaling(p, cropRange, leftTop, rightBottom, varargin) 
%
% Adapted from DJK_cropImages_3colors. However, if a fluorimage was not 2x2
% binned, it will be resized to 2x2 and the 'Binning' info will be changed.
% (Otherwise, the cropping cuts big parts of the image away).
% The program assumes phase contrast images to be 1x1 binned!

% DJK_cropImages creates crops of all images. Coordinates given must be even (in
% order for the fluor images to be cropped correctly)
%
% crops up to 3 fluorescence pictures as well
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
% Overwrite any schnitzcells parameters/defaults given optional fields/values
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
Dphase2   = dir([p.imageDir, [p.movieName,'-p-*2-*.tif'] ]);    %changed for Frenchstyle. (NW 21.1.2012. original: [p.movieName,'-p-*2-*.tif'] )
% All fluor1 images                                        % NW 11/12/02
Dfluor1=dir([p.imageDir, [p.movieName,sprintf('-%s-*.tif',p.fluor1)] ]);
% All fluor2 images                                        % NW 11/12/02
Dfluor2=dir([p.imageDir, [p.movieName,sprintf('-%s-*.tif',p.fluor2)] ]);
% All fluor3 images                                        % NW 11/12/02
Dfluor3=dir([p.imageDir, [p.movieName,sprintf('-%s-*.tif',p.fluor3)] ]);


% All fluor y images
%DfluorY = dir([p.imageDir, [p.movieName,'-y-*.tif'] ]);
% All fluor r images
%DfluorR = dir([p.imageDir, [p.movieName,'-r-*.tif'] ]);

if isempty(Dphase2)
    disp('You seem to have 1 phase image per frame');
    if ~isempty(Dfluor1)                                  % NW 11/12/02
        message=sprintf('You seem to have a fluor1 image (%s) for every %2.2f frames', p.fluor1,length(DphaseAll)/length(Dfluor1));
        disp(message);
    end
    if ~isempty(Dfluor2)
        message=sprintf('You seem to have a fluor2 image (%s) for every %2.2f frames', p.fluor2,length(DphaseAll)/length(Dfluor2));
        disp(message);
    end
    if ~isempty(Dfluor3)
        message=sprintf('You seem to have a fluor3 image (%s) for every %2.2f frames', p.fluor3,length(DphaseAll)/length(Dfluor3));
        disp(message);
    end
    
    %if ~isempty(DfluorY)
    %    disp(['You seem to a fluor image for every ' num2str(length(DphaseAll)/length(DfluorY)) ' frames']);
    %end 
    %if ~isempty(DfluorR)
    %    disp(['You seem to a fluor2 image for every ' num2str(length(DphaseAll)/length(DfluorY)) ' frames']);
    %end 
else 
    disp(['You seem to have ' num2str(length(DphaseAll)/length(Dphase2)) ' phase images per frame']);
    if ~isempty(Dfluor1)                                  % NW 11/12/02
        message=sprintf('You seem to have a fluor1 image (%s) for every %2.2f frames', p.fluor1,length(Dphase2)/length(Dfluor1));
        disp(message);
    end
    if ~isempty(Dfluor2)
        message=sprintf('You seem to have a fluor2 image (%s) for every %2.2f frames', p.fluor2,length(Dphase2)/length(Dfluor2));
        disp(message);
    end
    if ~isempty(Dfluor3)
        message=sprintf('You seem to have a fluor3 image (%s) for every %2.2f frames', p.fluor3,length(Dphase2)/length(Dfluor3));
        disp(message);
    end
    
    
    
    %if ~isempty(DfluorY)
    %    disp(['You seem to a fluor image for every ' num2str(length(Dphase2)/length(DfluorY)) ' frames ']);
    %end 
    %if ~isempty(DfluorY)
    %    disp(['You seem to a fluor2 image for every ' num2str(length(Dphase2)/length(DfluorY)) ' frames ']);
    %end 
end

if isempty(Dfluor1)                                  % NW 11/12/02
    disp('You seem to have no fluor1 images');
end
if isempty(Dfluor1)
    disp('You seem to have no fluor2 images');
end
if isempty(Dfluor3)
    disp('You seem to have no fluor3 images');
end

%if isempty(DfluorY)
%    disp('You seem to have no fluor images');
%end
%if isempty(DfluorR)
%    disp('You seem to have no fluor2 images');
%end
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
    %DfluorYRange = dir([p.imageDir, [p.movieName, '-y-*', str3(fr) ,'.tif'] ]);
    %DfluorRRange = dir([p.imageDir, [p.movieName, '-r-*', str3(fr) ,'.tif'] ]);
    Dfluor1Range = dir([p.imageDir, [p.movieName, sprintf('-%s-*',p.fluor1), str3(fr) ,'.tif'] ]);  % NW 11/12/02
    Dfluor2Range = dir([p.imageDir, [p.movieName, sprintf('-%s-*',p.fluor2), str3(fr) ,'.tif'] ]);  % NW 11/12/02
    Dfluor3Range = dir([p.imageDir, [p.movieName, sprintf('-%s-*',p.fluor3), str3(fr) ,'.tif'] ]);  % NW 11/12/02
    
   
    % crop phase images
    for i = [1:length(DphaseRange)], % go over each phase image of this frame
        % read image
        im_original = imread([p.imageDir DphaseRange(i).name]); 
        % get size of original image NW2012-10
        sizeph1=size(im_original,1);
        sizeph2=size(im_original,2);
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
    
    %crop fluor images (if existent)
    for i = [1:length(Dfluor1Range)],  % go over each fluor1 image of this frame (NW 11/12/02)
        % read image
        im_original = imread([p.imageDir Dfluor1Range(i).name]);
        % get size of fluor image NW2012-10
        sizefluor1=size(im_original,1);
        sizefluor2=size(im_original,2);
        % check if binning==2 and if not, rescale image
        scale1=sizeph1/sizefluor1;
        scale2=sizeph2/sizefluor2;
        if scale1~=scale2
            disp('Binning seems different for x and y axis. Can''t deal with that. Maybe try DJK_cropImages_3colors.')
        elseif scale1~=2
            im_original=imresize_old(im_original,scale1/2);
        end
        % get image info
        im_info = imfinfo([p.imageDir Dfluor1Range(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        % if image had to be rescaled, change the 'Binning' in
        % imageDescription NW2012-10
        if scale1~=2
            x=strfind(im_description,'Binning');
            descr(x+9:x+13)=['2 x 2'];
        end
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor1Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        if scale1==2
            disp(['Written: ' im_crop_filename]);
        else
            disp(['Rescaled to 2x2 Binning and Written: ' im_crop_filename]);
        end
        %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
    for i = [1:length(Dfluor2Range)],  % go over each fluor2 image of this frame (NW 11/12/02)
        % read image
        im_original = imread([p.imageDir Dfluor2Range(i).name]); 
        % get size of fluor image NW2012-10
        sizefluor1=size(im_original,1);
        sizefluor2=size(im_original,2);
        % check if binning==2 and if not, rescale image
        scale1=sizeph1/sizefluor1;
        scale2=sizeph2/sizefluor2;
        if scale1~=scale2
            disp('Binning seems different for x and y axis. Can''t deal with that. Maybe try DJK_cropImages_3colors.')
        elseif scale1~=2
            im_original=imresize_old(im_original,scale1/2);
        end
        % get image info
        im_info = imfinfo([p.imageDir Dfluor2Range(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        if scale1~=2
            x=strfind(im_description,'Binning');
            descr(x+9:x+13)=['2 x 2'];
        end
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor2Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
       if scale1==2
            disp(['Written: ' im_crop_filename]);
        else
            disp(['Rescaled to 2x2 Binning and Written: ' im_crop_filename]);
       end
       %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
    for i = [1:length(Dfluor3Range)],  % go over each fluor3 image of this frame (NW 11/12/02)
        % read image
        im_original = imread([p.imageDir Dfluor3Range(i).name]); 
        % get size of fluor image NW2012-10
        sizefluor1=size(im_original,1);
        sizefluor2=size(im_original,2);
        % check if binning==2 and if not, rescale image
        scale1=sizeph1/sizefluor1;
        scale2=sizeph2/sizefluor2;
        if scale1~=scale2
            disp('Binning seems different for x and y axis. Can''t deal with that. Maybe try DJK_cropImages_3colors.')
        elseif scale1~=2
            im_original=imresize_old(im_original,scale1/2);
        end
        % get image info
        im_info = imfinfo([p.imageDir Dfluor3Range(i).name]);
        % this image info will be added to crop
        im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
        if scale1~=2
            x=strfind(im_description,'Binning');
            descr(x+9:x+13)=['2 x 2'];
        end
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor3Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        if scale1==2
            disp(['Written: ' im_crop_filename]);
        else
            disp(['Rescaled to 2x2 Binning and Written: ' im_crop_filename]);
        end
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
