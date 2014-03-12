% DJK_analyzeFluorBackground 
%
% 
%
% OUTPUT
%
% REQUIRED ARGUMENTS:
%
% OPTIONAL ARGUMENTS:
% 'manualRange'   These frames will be treated
%
% 'DJK_saveDir'   Folder where tiff files are saved
%                 (default: '\analysis\fluor\')
%
% 'TIFFonly' = 1  does nothing with seg files, only creates a y2 TIFF 
%                 manualRange needs to be given, else might give
%                 error when loading seg file
%
% 'onScreen' = 1  will not automatically save and close images, but ask
%

function DJK_analyzeFluorBackground(p,varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_analyzeFluorBackground';

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
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]); return;
  end
end

% If explicit TIFFonly is not given, TIFFonly is false (0)
if ~existfield(p,'TIFFonly')
  p.TIFFonly = 0;
end

% Just in case it has not been set yet
if ~existfield(p,'onScreen')
  p.onScreen = 0;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Open file to write results to
%--------------------------------------------------------------------------
fid = fopen([p.DJK_saveDir p.movieName '-fluorBackground.txt'],'wt');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Let know what is happening
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SEG FILES AND GET BACKGROUND DATA
%--------------------------------------------------------------------------
all_frames = [];
all_yback = [];
all_yback2 = [];
all_yback2Alt = [];
all_yback3 = [];
all_ybackAlt = [];
all_ybackNew = [];
all_ybackNew2 = [];
all_ybackNew3 = [];


% loop over frames
count = 0;
for frameNum = p.manualRange
  %------------------------------------------------------------------------
  % LOAD DATA
  %------------------------------------------------------------------------
  % load complete fluor image
  yname= [p.imageDir,p.movieName,'-y-',str3(frameNum),'.tif'];
  if exist(yname)==2
    yimage = imread(yname);
  else
    dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-y-' str3(frameNum) '.tif in ' p.imageDir]);
    continue;
  end

  % load complete fluor image
  y2name= [p.analysisDir 'fluor\y2\' p.movieName '-y2-' str3(frameNum) '.tif'];
  if exist(y2name)==2
    y2image = imread(y2name);
  else
    dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-y2-' str3(frameNum) '.tif in ' p.analysisDir 'fluor\y2\']);
    continue;
  end

  if ~p.TIFFonly
    % load complete fluor image
    y3name= [p.analysisDir 'fluor\y3\' p.movieName '-y3-' str3(frameNum) '.tif'];
    if exist(y3name)==2
      y3image = imread(y3name);
    else
      dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-y3-' str3(frameNum) '.tif in ' p.analysisDir 'fluor\y3\']);
      continue;
    end

    % load segmentation file
    filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];
    clear yreg rect ybinning yback phaseFullSize gainy expty Lc LNsub;
    yreg = [];
    load(filename);
    dispAndWrite(fid, [' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);

    % load fluor file
    filename = [p.tracksDir p.movieName 'Fluor' str3(frameNum)];
    load(filename);
    dispAndWrite(fid, ['       -> loaded ' p.movieName 'Fluor' str3(frameNum) ' in ' p.tracksDir]);

    % get segmented part correct
    if ~exist('Lc')      
      dispAndWrite(fid, ['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
      Lc = LNsub;
    end
    segcrop = +(Lc > 0);
    seg = zeros(phaseCropSize);
    seg(rect(1):rect(3),rect(2):rect(4)) = segcrop;
  else
    % Get sizes correct
    phaseCropSize = phaseFullSize;
    seg = zeros(phaseCropSize);
  end

  % still need to resize yimage
  yimage = imresize(yimage,ybinning,'nearest');
  y2image = imresize(y2image,ybinning,'nearest');
  %------------------------------------------------------------------------
 
  
  %------------------------------------------------------------------------
  % COLLECT BACKGROUND DATA
  %------------------------------------------------------------------------
  all_frames = [all_frames frameNum];
  all_yback = [all_yback yback];
  all_yback2 = [all_yback2 yback2];
  all_yback2Alt = [all_yback2Alt yback2Alt];
  all_yback3 = [all_yback3 yback3];
  all_ybackAlt = [all_ybackAlt ybackAlt];
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % GET NEW BACKGROUND DATA 1
  %------------------------------------------------------------------------
  seg_dilate = imdilate(seg,strel('disk',30));
  [rows, cols] = find(seg_dilate>0);
  row_min = min(rows);
  row_max = max(rows);
  col_min = min(cols);
  col_max = max(cols);
  rect_dilate = [row_min col_min row_max col_max];

  background1 = zeros(phaseCropSize);
  background1(row_min:row_max,col_min:col_max) = 1;
  background1(find(seg_dilate>0)) = 0;
  idx = find(background1>0);
    
  ybackNew = median(y3image(idx));
  all_ybackNew = [all_ybackNew ybackNew];
  %------------------------------------------------------------------------

  
  %------------------------------------------------------------------------
  % GET NEW BACKGROUND DATA 2 & 3
  %------------------------------------------------------------------------
  % first get boundary of microcolonie
  seg_filled = imfill(imclose(seg, strel('disk',1)));
  seg_filled_perim = bwperim(seg_filled);
  [perim_y, perim_x] = find(seg_filled_perim>0);
  
  % get min distance for each pixel away from seg
  [nonseg_y, nonseg_x] = find(seg_filled==0);
  nonseg_idx = find(seg_filled==0);
  
  %disp('starting calc optimized'); t=cputime; 
  for i = 1:length(perim_y)
    dist_squared = (nonseg_y-perim_y(i)).^2 + (nonseg_x-perim_x(i)).^2;
    if i==1, min_dist_squared = dist_squared;
    else, min_dist_squared = min(min_dist_squared, dist_squared);
    end
  end
  %disp(['finished calc optimized in ' num2str(cputime-t)]); pause(0.2);
  % keyboard;

  pixels_squared_per_micro = (1/p.micronsPerPixel)^2;
  
  % 0 tot 1 um etc
  for i = 1:7
    pixels2 = y3image( nonseg_idx(find(min_dist_squared>=(i-1)*pixels_squared_per_micro & min_dist_squared<=i*pixels_squared_per_micro)) );
    pixels3 = y2image( nonseg_idx(find(min_dist_squared>=(i-1)*pixels_squared_per_micro & min_dist_squared<=i*pixels_squared_per_micro)) );
    if length(pixels2) > 0, ybackNew2(i) = median(pixels2); ybackNew3(i) = median(pixels3);
    else,                   ybackNew2(i) = 0; ybackNew3(i) = 0;
    end
  end
  count = count + 1;
  all_ybackNew2(:,count) = [ybackNew2(1); ybackNew2(2); ybackNew2(3); ybackNew2(4); ybackNew2(5); ybackNew2(6); ybackNew2(7); ];
  all_ybackNew3(:,count) = [ybackNew3(1); ybackNew3(2); ybackNew3(3); ybackNew3(4); ybackNew3(5); ybackNew3(6); ybackNew3(7); ];
  %------------------------------------------------------------------------
end 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PRINT BACKGROUND DATA
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['frame yback yback2 yback2Alt yback3 ybackAlt ybackNew 0-1 1-2 2-3 3-4 4-5 5-6 6-7']);
for i = 1:length(all_frames)
  dispAndWrite(fid, [ str3(all_frames(i)) ' ' ...
                      str3(all_yback(i)) ' ' ...
                      str3(all_yback2(i)) ' ' ...
                      str3(all_yback2Alt(i)) ' ' ...
                      str3(all_yback3(i)) ' ' ...
                      str3(all_ybackAlt(i)) ' ' ...
                      str3(all_ybackNew(i)) ' ' ...
                      str3(all_ybackNew2(1,i)) ' ' ...
                      str3(all_ybackNew2(2,i)) ' ' ...
                      str3(all_ybackNew2(3,i)) ' ' ...
                      str3(all_ybackNew2(4,i)) ' ' ...
                      str3(all_ybackNew2(5,i)) ' ' ...
                      str3(all_ybackNew2(6,i)) ' ' ...
                      str3(all_ybackNew2(7,i)) ' ' ...
                      str3(all_ybackNew3(1,i)) ' ' ...
                      str3(all_ybackNew3(2,i)) ' ' ...
                      str3(all_ybackNew3(3,i)) ' ' ...
                      str3(all_ybackNew3(4,i)) ' ' ...
                      str3(all_ybackNew3(5,i)) ' ' ...
                      str3(all_ybackNew3(6,i)) ' ' ...
                      str3(all_ybackNew3(7,i)) ' ' ]);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% NORMALIZE FLUOR DATA SAME WAY AS FOR SCHNITZES
%--------------------------------------------------------------------------
all_yback = normalize( all_yback, 0, expty, ybinning, p.micronsPerPixel, 0);
all_yback2 = normalize( all_yback2, 0, expty, ybinning, p.micronsPerPixel, 0);
all_yback2Alt = normalize( all_yback2Alt, 0, expty, ybinning, p.micronsPerPixel, 0);
all_yback3 = normalize( all_yback3, 0, expty, ybinning, p.micronsPerPixel, 0);
all_ybackAlt = normalize( all_ybackAlt, 0, expty, ybinning, p.micronsPerPixel, 0);
all_ybackNew = normalize( all_ybackNew, 0, expty, ybinning, p.micronsPerPixel, 0);
all_ybackNew2 = normalize( all_ybackNew2, 0, expty, ybinning, p.micronsPerPixel, 0);
all_ybackNew3 = normalize( all_ybackNew3, 0, expty, ybinning, p.micronsPerPixel, 0);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Show Plot
%--------------------------------------------------------------------------
% make figure
scrsz = get(0, 'ScreenSize');
fig111 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
hold on;

subplot(2,1,1);
plot(all_frames, all_yback, 'b-', all_frames, all_ybackAlt,'k:'); 
legend('yback','ybackAlt',-1);
ylim([10 15]);

subplot(2,1,2);
plot(all_frames, all_yback2, 'b-', all_frames, all_yback2Alt,'k:', all_frames, all_yback3,'b--',all_frames, all_ybackNew,'r-.'); 
legend('yback2','yback2Alt','yback3','ybackNew',-1);
ylim([0 5]);
xlabel('Frames');

fig112 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
hold on;

subplot(2,1,1);
plot( all_frames, all_ybackNew2(1,:), 'g-', ...
      all_frames, all_ybackNew2(2,:), 'r-', ...
      all_frames, all_ybackNew2(3,:), 'b-', ...
      all_frames, all_ybackNew2(4,:), 'g--', ...
      all_frames, all_ybackNew2(5,:), 'r--', ...
      all_frames, all_ybackNew2(6,:), 'b--', ...
      all_frames, all_ybackNew2(7,:), 'g-.'); 
legend('ybackNew2 0 to 1 um','ybackNew2 1 to 2 um','ybackNew2 2 to 3 um','ybackNew2 3 to 4 um','ybackNew2 4 to 5 um','ybackNew2 5 to 6 um','ybackNew2 6 to 7 um',-1);
ylim([0 5]);

subplot(2,1,2);
plot( all_frames, all_ybackNew3(1,:), 'g-', ...
      all_frames, all_ybackNew3(2,:), 'r-', ...
      all_frames, all_ybackNew3(3,:), 'b-', ...
      all_frames, all_ybackNew3(4,:), 'g--', ...
      all_frames, all_ybackNew3(5,:), 'r--', ...
      all_frames, all_ybackNew3(6,:), 'b--', ...
      all_frames, all_ybackNew3(7,:), 'g-.'); 
legend('ybackNew3 0 to 1 um','ybackNew3 1 to 2 um','ybackNew3 2 to 3 um','ybackNew3 3 to 4 um','ybackNew3 4 to 5 um','ybackNew3 5 to 6 um','ybackNew3 6 to 7 um',-1);
ylim([0 5]);
xlabel('Frames');

% Ask to save the figure
if p.onScreen
  saveFigInput = questdlg('Save Figures?','Save Figures?','Yes','Yes and Close','No','Yes');
  pause(0.2);
else
  saveFigInput='Yes and Close';
end

if (upper(saveFigInput(1))=='Y')
  saveas(fig111,[p.DJK_saveDir p.movieName '-fluorBackground1.fig']);
  saveSameSize(fig111,'file',[p.DJK_saveDir p.movieName '-fluorBackground1.png'], 'format', 'png');
  saveas(fig112,[p.DJK_saveDir p.movieName '-fluorBackground2.fig']);
  saveSameSize(fig112,'file',[p.DJK_saveDir p.movieName '-fluorBackground2.png'], 'format', 'png');
  if (strcmp(saveFigInput,'Yes and Close'))
    close(fig111); close(fig112);
    pause(0.2);
  end
  dispAndWrite(fid, ['-------------------------------------------------']);
  dispAndWrite(fid, [' * Saved plot backgrund fluor over time in ' p.movieName '-fluorBackground1.png']);
  dispAndWrite(fid, [' * Saved plot backgrund fluor over time in ' p.movieName '-fluorBackground2.png']);
end
%--------------------------------------------------------------------------
  

%--------------------------------------------------------------------------
% Close file to write results to
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
if fid~=0
  fclose(fid);
end
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize performs new normalization
function yregNorm = normalize( yreg, ...            % gfp image in camera units
                               yback, ...           % background in camera units
                               expty, ...           % exposure time in ms
                               ybinning, ...        % camera bin size (2)
                               micronsPerPixel, ... % micron per pixel (0.04065)
                               autoFluor)           % autofluorescence of cells
if (expty=='empty')
  disp('Warning! no exposure time specified for fluor image');
  expty = 1000; 
end
yregNorm=double(yreg)-double(yback);
yregNorm=yregNorm/(expty * ybinning^2 * micronsPerPixel^2);
yregNorm=yregNorm-autoFluor; % in units of molecules?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%