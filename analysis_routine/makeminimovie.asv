function newp = makeminimovie(p,varargin) 
% MAKEMINIMOVIE creates a sub-movie within the given movie around clicked points
% 
%   MAKEMINIMOVIE creates a new movie that consists of mini sub-images 
%   extracted from frames of a larger movie.  The routine allows the user to 
%   click to define the center of each new mini movie frame, by showing the 
%   user successive frames of the given movie from the end of the movie to 
%   the start.  A new movie directory structure is created, consisting of 
%   sub-images of winSize pixels that are written to the new mini-movie's 
%   image directory.  MAKEMINIMOVIE returns a new schnitzcells movie parameter 
%   structure so the new mini movie can be analyzed.
%
% Optional Property/Value inputs:
%
%   'miniMovieName'   New name for movie, by default [p.movieName+'-mini-01']; 
%                     this becomes the name of the new movie, newp.movieName.
% 
%   'winSize'         Scalar dimension of sub-images extracted around clicks.
%                     Default is 200 (resulting in 200x200 sub-images).
% 
%   'miniRange'       Defines a specific range of frame numbers to extract. 
%                     By default, new movie will contain as many frames as 
%                     the given original movie.
% 
% NOTE: Still need to add code to handle movies with 2 or more phase slices 
% better. This routine should rewrite ALL phase images, but presently only 
% rewrites one, and sets numphaseslices of new movie to 1.


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
   (mod(nargin,2) == 0) | ...
   (~isSchnitzParamStruct(p))
  errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> makeminimovie:',...
      '    Invalid input arguments.',...
      '    Try "help makeminimovie".');
  error(errorMessage);
end

%-------------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values 
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
          'Error using ==> absaugen:',...
          '    Invalid property ', num2str(varargin{i}), ...
          ' is not (needs to be) a string.',...
          '    Try "help absaugen".');
      error(errorMessage);
    end
    fieldName = schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end


if ~existfield(p,'miniMovieName')
  p.miniMovieName = [p.movieName,'-mini-01'];
end

if ~existfield(p,'winSize')
  p.winSize = 200;
end

% If explicit frame range not given, figure it out by determining which frames 
% exist (based on image file names)
imprefix = [p.movieName '-y-'];
D = dir([p.imageDir, imprefix, '*.tif']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.tif')-3;
if ~existfield(p,'miniRange')
  D = dir([p.imageDir, imprefix, '*.tif']);
  [S,I] = sort({D.name}');
  D = D(I);
  frameNumStrings = char(S);
  p.miniRange = str2num(frameNumStrings(:,numpos:numpos+2))';
end


% Internally, code is set up to use the mini frame range in descending order
sortedMiniRange = sort(p.miniRange); % in case they give it to us descending!
descendingMiniRange = sortedMiniRange(end:-1:1);
i = length(p.miniRange);  % i = row index into frame info, NOT frame number!
p.numphaseslices = 1; %%%%%%%%%%%ADD SJT
for frame=descendingMiniRange
    if p.numphaseslices == 1
        imp = imread([p.imageDir,p.movieName,'-p-',str3(frame),'.tif']);
    else
        % assumes phase image 2 is the one to extract from!
        imp = imread([p.imageDir,p.movieName,'-p-2-',str3(frame),'.tif']);
    end
    imc = imread([p.imageDir,p.movieName,'-y-',str3(frame),'.tif']);
    imp = imresize(imp,0.5);
    imshow(makergb(imp,imc));
    frame,
    w = waitforbuttonpress;
    if w == 0,
    
        F(i) = frame;
        xy = 2*round(get(gca,'CurrentPoint'));
        X(i) = xy(1,1);
        Y(i) = xy(1,2);
        [frame X(i) Y(i)]
        i = i - 1;

    end;
end;

clickdata = [F' X' Y'];


frames = clickdata(:,1);
X = clickdata(:,2);
Y = clickdata(:,3);

% X = X*2;
% Y = Y*2;

% create a new movie with only one phase slice
newp = initschnitz(p.miniMovieName,p.movieDate,p.movieKind,'rootDir',p.rootDir,'numphaseslices',1);

if p.numphaseslices == 1
    basepath = [p.imageDir,p.movieName,'-p-'];
else
    % assumes phase image 2 is the one to extract from!
    basepath = [p.imageDir,p.movieName,'-p-2-'];
end
outpath = [newp.imageDir,newp.movieName,'-p-'];
for i = 1:size(clickdata,1),
    
    readname = [basepath,str3(frames(i)),'.tif']
    
    imp = imread(readname);
    sz = size(imp);

    startx = round(max(1,X(i)-p.winSize/2));
    stopx = round(min(sz(2),X(i)+p.winSize/2));

    starty = round(max(1,Y(i)-p.winSize/2));
    stopy = round(min(sz(1),Y(i)+p.winSize/2));

    pwin = imp(starty:stopy,startx:stopx);
    
    imwrite(pwin,[outpath,str3(frames(i)),'.tif'],'tif');
    
end;
