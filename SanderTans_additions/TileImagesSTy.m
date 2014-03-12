function p = TileImagesST(p,varargin) 

% p = TileImagesSTy(p,'numphaseslices',1,'NrTilesX',5,'NrTilesY',5,'EnlargeFact',1.3,'tileRange',1:25);

% p = TileImagesSTy(p,'numphaseslices',1,'NrTilesX',7,'NrTilesY',5,'EnlargeFact',1.3,'tileRange',1:35);
%
% when printing, first resize screen to resemble A4, do 'print review',
% 'page set up','fix aspext ratio', en 'fill page'
%
%


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

global ourfig origfig reducefig

tilefig=2;
reducefig=3;

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

F=[];X=[];Y=[];
dx=0.09;dy=0.06;
ddx=dx+0.05;ddy=dy+0.01;
nrTilesX=5;
nrTilesY=7;

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

if ~existfield(p,'numphaseslices')
  p.numphaseslices = 1;
end

if ~existfield(p,'NrTilesX')
  p.NrTilesX = 5;
end

if ~existfield(p,'NrTilesY')
  p.NrTilesX = 5;
end

if ~existfield(p,'EnlargeFact')
  p.EnlargeFact = 1.2;
end

% If explicit frame range not given, figure it out by determining which frames 
% exist (based on image file names)
imprefix = [p.movieName '-y-'];
D = dir([p.imageDir, imprefix, '*.tif']);
[S,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.tif')-3;
if ~existfield(p,'tileRange')
  D = dir([p.imageDir, imprefix, '*.tif']);
  [S,I] = sort({D.name}');
  D = D(I);
  frameNumStrings = char(S);
  p.miniRange = str2num(frameNumStrings(:,numpos:numpos+2))';
end

% bgDir=[p.rootDir,'Background\bg2.tif'];
% BG=imread(bgDir);

% Internally, code is set up to use the mini frame range in descending order
sortedMiniRange = sort(p.tileRange); % in case they give it to us descending!
descendingMiniRange = sortedMiniRange(end:-1:1);
descendingMiniRange = p.tileRange; % ADD: do not sort SJT
descendingMiniRange;
p.tileRange;
i=1;
imax = length(descendingMiniRange);  % i = row index into frame info, NOT frame number!
% p.numphaseslices = 2; %%%%%%%%%%%ADD SJT
% for frame=descendingMiniRange
figure(tilefig);
% set(gca,'ActivePositionProperty','OuterPosition');
while i<=imax
    frame=descendingMiniRange(i);
    if p.numphaseslices == 1
        imp = imread([p.imageDir,p.movieName,'-p-',str3(frame),'.tif']);
    else
        % assumes phase image 2 is the one to extract from!
        imp = imread([p.imageDir,p.movieName,'-p-2-',str3(frame),'.tif']);
    end
    imy = imread([p.imageDir,p.movieName,'-y-',str3(frame),'.tif']);
    %     imy = imy + 600 - BG;
    %    imp = imresize(imp,0.5);
        %iptsetpref('ImshowAxesVisible','on','ImshowBorder','loose' );
    subplot(p.NrTilesX,p.NrTilesY,i);
%     x=ddx+ddx*fix((i-1)/5);
%     y=1-1*ddy-(ddy*rem(i-1,5));

%     subplot('Position',[y x dy dx]);
     axis off;
%     [i get(gca,'OuterPosition')]
%     
    imshow(makergb(imresize(imp,0.5),imy));
    d=get(gca,'Position');
    c=p.EnlargeFact;
      set(gca,'Position',[d(1) d(2) c*d(3) c*d(4)]);    
%      set(gca,'Position',(get(gca,'Position')+get(gca,'TightInset')));
       %iptsetpref('ImshowAxesVisible','on');
     %    iptsetpref('ImshowBorder','loose');
    %     frame
    i=i+1;
end

TitleString1=['File:  ' p.movieDir '  RANGE ' num2str(min(p.tileRange)) '   to   ' num2str(max(p.tileRange))];
TagTitle = axes('Position',[0.0 -0.0 0.8 0.04]);
axis off;
text(0.2,0.7,[TitleString1],'FontSize',8);

