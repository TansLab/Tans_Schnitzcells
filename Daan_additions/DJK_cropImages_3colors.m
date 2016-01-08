function DJK_cropImages_3colors(p, cropRange, leftTop, rightBottom, varargin) 
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
% Optional inputs
%   p.resizePhase     This is to handle very exceptional cases. If you have
%                     a phaseimage with an incorrect size (e.g. because you
%                     have binned your phase image to another binning than
%                     the required 1x1 binning), you can correct that by
%                     setting this parameter. E.g., if p.resizePhase = 2
%                     than the cropped phase image will become twice as
%                     large. 

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings

numRequiredArgs = 4;
functionName = 'DJK_cropImages';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
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
% overwrite any schnitzcells parameters/defaults given optional fields/values
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
if isempty(Dfluor2)
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
    	
    DphaseRange = dir([p.imageDir, [p.movieName, '-p-*-', str3(fr) ,'.tif'] ]); % fix to handle frame nrs >999 MW 2014/6/3

    %DfluorYRange = dir([p.imageDir, [p.movieName, '-y-*', str3(fr) ,'.tif'] ]);
    %DfluorRRange = dir([p.imageDir, [p.movieName, '-r-*', str3(fr) ,'.tif'] ]);
    Dfluor1Range = dir([p.imageDir, [p.movieName, sprintf('-%s-',p.fluor1), str3(fr) ,'.tif'] ]);  % NW 11/12/02 + MW 2014/6, fix to handle nrs >999
    Dfluor2Range = dir([p.imageDir, [p.movieName, sprintf('-%s-',p.fluor2), str3(fr) ,'.tif'] ]);  % NW 11/12/02 + MW 2014/6
    Dfluor3Range = dir([p.imageDir, [p.movieName, sprintf('-%s-',p.fluor3), str3(fr) ,'.tif'] ]);  % NW 11/12/02 + MW 2014/6
    
    if isempty(DphaseRange)
        disp(['***' 10 'Searching for img with : ''', p.movieName, '-p-*-', str3(fr) ,'.tif', '''']);
        disp(['Looking in directory ' p.imageDir]);
        error(['Couldn''t find image, perhaps img dir empty/frame missing?' 10 '***'])        
    end
   
    % crop phase images
    for i = [1:length(DphaseRange)], % go over each phase image of this frame

        % read image       
        im_original = imread([p.imageDir DphaseRange(i).name]); 


        % Resize the image. This is an option to be used for exceptional
        % cases. E.g. when the phase image has been accidentally binned,
        % other parts of the code will not function poroperly. Here, the
        % phase image can be converted to a resolution that would
        % correspond to 1x1 binning. (Currently, a resolution of
        % 2048x2048.)
        if isfield(p, 'resizePhase')
            im_original = imresize(im_original, p.resizePhase);
        end        
        
        % crop the image
        im_crop = im_original( leftTop(2):rightBottom(2), leftTop(1):rightBottom(1));
       
        % read iminfo
        im_info = imfinfo([p.imageDir DphaseRange(i).name]);
       
        % this image info will be added to crop        
        if strcmp(p.softwarePackage, 'micromanager')
            im_description=DE_adjustiminfo(p, DphaseRange(i).name);
        elseif strcmp(p.softwarePackage, 'metamorph')
            % old version only allowed to add one extra field.
           % im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];

            % NW 2014-11. Add Date&Software info which are otherwise
            % lost. Some weird matlab behaviour requires the extra
            % fields to be converted to strings ('char' not enough!)
            % first. Otherweise only the first field (Date) is stored
            % in image info - yes even though the im_description
            % contains both fields only the first field will be stored!
            im_description = [im_info.ImageDescription];
            extraInfo1=sprintf(['DateTime: ' im_info.DateTime]);
            extraInfo2=sprintf([' Software: ' im_info.Software]);
            im_description=[im_description extraInfo1 extraInfo2];
        else
            error('No software package was set in p.softwarePackage or not recognized.')
        end

        
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(DphaseRange(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written frame ' num2str(fr) ' : ' im_crop_filename]);
        %disp([ num2str(leftTop(2)) ':' num2str(rightBottom(2)) '-' num2str(leftTop(1)) ':' num2str(rightBottom(1))]);
    end
        
    %crop fluor images (if existent)
    for i = [1:length(Dfluor1Range)],  % go over each fluor1 image of this frame (NW 11/12/02)
        
        % read image
        im_original = imread([p.imageDir Dfluor1Range(i).name]);
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % read iminfo
        im_info = imfinfo([p.imageDir Dfluor1Range(i).name]);

        
         % this image info will be added to crop
        if strcmp(p.softwarePackage,'micromanager')
            im_description=DE_adjustiminfo(p, Dfluor1Range(i).name);
        elseif strcmp(p.softwarePackage,'metamorph')
               % Add Date&Software info. Details see above.
               im_description = [im_info.ImageDescription];
               extraInfo1=sprintf(['DateTime: ' im_info.DateTime]);
               extraInfo2=sprintf([' Software: ' im_info.Software]);
               im_description=[im_description extraInfo1 extraInfo2];
        else
            error('No software package was set in p.softwarePackage or not recognized.')
        end        
    
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor1Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written frame ' num2str(fr) ' : ' im_crop_filename]);
        %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
    
    for i = [1:length(Dfluor2Range)],  % go over each fluor2 image of this frame (NW 11/12/02)
        
        % read image
        im_original = imread([p.imageDir Dfluor2Range(i).name]); 
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
        % read iminfo
        im_info = imfinfo([p.imageDir Dfluor2Range(i).name]);
        
          % this image info will be added to crop
        if strcmp(p.softwarePackage,'micromanager')
            im_description=DE_adjustiminfo(p, Dfluor2Range(i).name);
        elseif strcmp(p.softwarePackage,'metamorph')
                % Add Date&Software info. Details see above.
               im_description = [im_info.ImageDescription];
               extraInfo1=sprintf(['DateTime: ' im_info.DateTime]);
                extraInfo2=sprintf([' Software: ' im_info.Software]);
               im_description=[im_description extraInfo1 extraInfo2];
        else
            error('No software package was set in p.softwarePackage or not recognized.')
        end
        
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor2Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written frame ' num2str(fr) ' : ' im_crop_filename]);
        %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
    end
    
    for i = [1:length(Dfluor3Range)],  % go over each fluor3 image of this frame (NW 11/12/02)
        % read image
        im_original = imread([p.imageDir Dfluor3Range(i).name]); 
        % crop the image
        im_crop = im_original( ((leftTop(2)+1)/2):(rightBottom(2)/2), ((leftTop(1)+1)/2):(rightBottom(1)/2));
         % read iminfo
        im_info = imfinfo([p.imageDir Dfluor3Range(i).name]);
        
         % this image info will be added to crop
        if strcmp(p.softwarePackage,'micromanager')
                im_description=DE_adjustiminfo(p, Dfluor3Range(i).name);
        elseif strcmp(p.softwarePackage,'metamorph')
               % Add Date&Software info. Details see above.
               im_description = [im_info.ImageDescription];
               extraInfo1=sprintf(['DateTime: ' im_info.DateTime]);
                extraInfo2=sprintf([' Software: ' im_info.Software]);
               im_description=[im_description extraInfo1 extraInfo2];
        else
            error('No software package was set in p.softwarePackage or not recognized.')
        end

        
        % write image data
        im_crop_filename = [cropDirImageDir regexprep(Dfluor3Range(i).name, p.movieName, p.cropName)];
        imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
        disp(['Written frame ' num2str(fr) ' : ' im_crop_filename]);
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
