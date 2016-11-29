function p = DJK_analyzeSeg(p, varargin);
% DJK_analyzeSeg analyses the segmentation files.
%
% Reports on:
%  * whether there are frames whose segmentation is not approved
%  * drops in cell number in subsequent frames
%  * saves number of cells in each frame to text file
%  * plots nrCells vs Time
%  * plots length, width and area of cells
%
%  Will save plots and close them, unless onScreen=1 is used.
%
% OPTIONAL ARGUMENTS:
%   'manualRange' allows to analyze a subset of frames (standard: all framse)
%
%   'fitRange' allows to fit a data of a subset of frames (standard: manualRange)
%
%   'micronsPerPixel' sets the pixel size (0.04065 standard)
%
%   'onScreen' = 1 will not automatically save and close images, but ask
%
%   'DJK_saveDir'       Directory where images will be saved. Defaults to
%                       "p.analysisDir 'segmentation\' 'fitRange1_200\'"
%

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1;
functionName = 'DJK_analyzeSeg';

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

% If explicit fitRange is not given, use manualRange
if ~existfield(p,'fitRange')
  p.fitRange = p.manualRange;
end

% Just in case it has not been set yet
if ~existfield(p,'micronsPerPixel')
  p.micronsPerPixel = 0.04065;
end

% Just in case it has not been set yet
if ~existfield(p,'onScreen')
  p.onScreen = 0;
end

% Just in case it has not been set yet
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'segmentation' filesep 'fitRange' num2str(p.fitRange(1)) '_' num2str(p.fitRange(end)) filesep];
end

% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Open file to write results to
%--------------------------------------------------------------------------
fid = fopen([p.DJK_saveDir p.movieName '-segmentation.txt'],'wt');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOAD SEGMENTATION FILES AND TAKE OUT DATA TO CHECK
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
dispAndWrite(fid, [' * fitRange is ' num2str(length(p.fitRange)) ' frames from ' num2str(p.fitRange(1)) ' to ' num2str(p.fitRange(end))]);
dispAndWrite(fid, [' * Saving files in : ' p.DJK_saveDir]);
dispAndWrite(fid, [' * Loading segmentation files...']);

i = 0;
% loop over frames
for frameNum = p.manualRange
  i = i+1;
  clear LNsub Lc timestamp;
  load([p.segmentationDir, p.movieName, 'seg', str3(frameNum)]);
 
  %  remember frames whose segmentation is not approved
  if ~(exist('Lc') == 1)
    unApproved(i) = 1;
    % will use Lc(corrected), so if not approved copy unapproved to Lc
    Lc = LNsub;
  else
    unApproved(i) = 0;
  end

  % Get number of segmented cells in frame
  NrCellsInFrame(frameNum) = max2(Lc);

  % Get time of frame  
  TimeOfFrame(frameNum) = timestamp;
  
  % Get length, width and area of cells in frame
  rp = regionprops(Lc,'majoraxislength','MinorAxisLength','Area');

  % calculate sums and averages
  sum_length = 0; sum_width = 0; sum_area = 0;
  for cellnum = 1:length(rp)
    sum_length = sum_length + rp(cellnum).MajorAxisLength;
    sum_width = sum_width + rp(cellnum).MinorAxisLength;
    sum_area = sum_area + rp(cellnum).Area;
  end
  length_sum(frameNum) = sum_length*p.micronsPerPixel;
  length_average(frameNum) = sum_length*p.micronsPerPixel/length(rp);
  width_sum(frameNum) = sum_width*p.micronsPerPixel;
  width_average(frameNum) = sum_width*p.micronsPerPixel/length(rp);
  area_sum(frameNum) = sum_area*p.micronsPerPixel*p.micronsPerPixel;
  area_average(frameNum) = sum_area*p.micronsPerPixel*p.micronsPerPixel/length(rp);
end;

% convert timestamp to minutes
TimeOfFrame = convert_2_minutes(TimeOfFrame);

% select data for fitting with fitRange
NrCellsInFrameFitting = NrCellsInFrame(p.fitRange);
TimeOfFrameFitting = TimeOfFrame(p.fitRange);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write cell number for each frame to file
%--------------------------------------------------------------------------
fid2 = fopen([p.DJK_saveDir p.movieName '-nrCells.txt'],'wt');
fprintf(fid2, '%s\n', ['frameNum' ' , ' 'nrCells']); 
for frameNum = p.manualRange
  fprintf(fid2, '%s\n', [num2str(frameNum) ' , ' num2str(NrCellsInFrame(frameNum))]);
end
fclose(fid2);
dispAndWrite(fid, [' * Saved Number of cells in each frame in ' p.movieName '-nrCells.txt']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Plots nrCells vs Time logarithmically
%--------------------------------------------------------------------------
% make plot
fig1 = figure;
semilogy(TimeOfFrame, NrCellsInFrame, 'ok');
hold on;
ylabel('# of cells');
xlabel('Time (mins)');
title([p.movieDate ' ' p.movieName ' nrCells-vs-time']);

% add fitted line
[muNrCells A0] = DJK_ExponentialFit(TimeOfFrameFitting/60, NrCellsInFrameFitting);
NrCells_fitted = A0*power(2,muNrCells/60*TimeOfFrameFitting);
semilogy(TimeOfFrameFitting, NrCells_fitted, '-','LineWidth',3);
text(0.05,0.9,['mu is ', num2str(muNrCells)],'sc');

% MW addition larger font size
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','normal')
set(gca,'FontSize',15)

% Ask to save the figure
if p.onScreen
  saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
  pause(0.2);
else
  saveFigInput='Yes and Close';
end
if (upper(saveFigInput(1))=='Y')
  saveas(fig1,[p.DJK_saveDir p.movieName '-nrCells-vs-time.fig']);
  saveSameSize(fig1,'file',[p.DJK_saveDir p.movieName '-nrCells-vs-time.eps'], 'format', 'epsc');
  saveas(fig1,[p.DJK_saveDir p.movieName '-nrCells-vs-time.tif'], 'tif');
  if (strcmp(saveFigInput,'Yes and Close'))
    close(fig1);
    pause(0.2);
  end
  dispAndWrite(fid, [' * Saved plot of number of cells vs time written in ' p.movieName '-nrCells-vs-time.png']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Plots nrCells vs Time logarithmically
%--------------------------------------------------------------------------
% make figure
screensize = get(0, 'ScreenSize');
fig2 = figure('Position', screensize);
hold on;

% make plot 1
subplot(3,2,1);
semilogy(TimeOfFrame, length_sum, 'ok');
ylabel('sum of cell lengths (um)');
xlabel('Time (mins)');
title([p.movieDate ' ' p.movieName ' CellSize-vs-time']);
% add fitted line
length_sumFitting = length_sum(p.fitRange);
[muLength A0] = DJK_ExponentialFit(TimeOfFrameFitting/60, length_sumFitting);
length_fitted = A0*power(2,muLength/60*TimeOfFrameFitting);
hold on;
semilogy(TimeOfFrameFitting, length_fitted, '-','LineWidth',3);
text(0.05,0.9,['mu is ', num2str(muLength)],'sc');

% make plot 2
subplot(3,2,2);
plot(TimeOfFrame, length_average, 'ok');
ylabel('average cell length (um)');
xlabel('Time (mins)');

% make plot 3
subplot(3,2,3);
semilogy(TimeOfFrame, area_sum, 'ok');
ylabel('sum of cell area (um^2)');
xlabel('Time (mins)');
% add fitted line
area_sumFitting = area_sum(p.fitRange);
[muVol A0] = DJK_ExponentialFit(TimeOfFrameFitting/60, area_sumFitting);
area_fitted = A0*power(2,muVol/60*TimeOfFrameFitting);
hold on;
semilogy(TimeOfFrameFitting, area_fitted, '-','LineWidth',3);
text(0.05,0.9,['mu is ', num2str(muVol)],'sc');

% make plot 4
subplot(3,2,4);
plot(TimeOfFrame, area_average, 'ok');
ylabel('average cell area (um^2)');
xlabel('Time (mins)');

% make plot 5
subplot(3,2,5);
semilogy(TimeOfFrame, width_sum, 'ok');
ylabel('sum of cell widths (um)');
xlabel('Time (mins)');

% make plot 6
subplot(3,2,6);
plot(TimeOfFrame, width_average, 'ok');
ylabel('average cell width (um)');
xlabel('Time (mins)');

% Ask to save the figure
if p.onScreen
  saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
  pause(0.2);
else
  saveFigInput='Yes and Close';
end
if (upper(saveFigInput(1))=='Y')
    saveas(fig2,[p.DJK_saveDir p.movieName '-CellSize-vs-time.fig']);
    saveSameSize(fig2,'file',[p.DJK_saveDir p.movieName '-CellSize-vs-time.eps'], 'format', 'epsc');
    saveas(fig2,[p.DJK_saveDir p.movieName '-CellSize-vs-time.tif'], 'tif');
    if (strcmp(saveFigInput,'Yes and Close'))
        close(fig2);
        pause(0.2);
    end
  dispAndWrite(fid, [' * Saved plot of cell size vs time written in ' p.movieName '-nrCells-vs-time.png']);
end
dispAndWrite(fid, ['-------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Show frames whose segmentation are not approved
%--------------------------------------------------------------------------
if (sum(unApproved)<1)
  dispAndWrite(fid, ['All frames are approved']);
else
  dispAndWrite(fid, [num2str(sum(unApproved)) ' frames are not approved:']);
  unapproved_str = []; % for printing
  for i=[1:length(unApproved)]
    if unApproved(i)
        unapproved_str = [unapproved_str ' ' str3(p.manualRange(i))];      
    end
  end
  dispAndWrite(fid,['* ' unapproved_str '.']);
end  
dispAndWrite(fid, ['-------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Check for drop in cell number in subsequent frames for approved frames
%--------------------------------------------------------------------------
dispAndWrite(fid, ['Checking for drops in cell number in subsequent approved frames:']);
nrCellsInPreviousFrame = 0;
for i=[1:length(unApproved)]
  frameNum = p.manualRange(i);
  if ~unApproved(i)
    if ( NrCellsInFrame(frameNum) < nrCellsInPreviousFrame )
      dispAndWrite(fid, [' * frame ', str3(p.manualRange(i-1)), ' -> ', str3(frameNum), ' : from ', ...
                        num2str(nrCellsInPreviousFrame), ' to ', num2str(NrCellsInFrame(frameNum)), ' cells']);
    end
  end
  nrCellsInPreviousFrame = NrCellsInFrame(frameNum);
end
dispAndWrite(fid, ['-------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Display some additional growth info
%--------------------------------------------------------------------------
dispAndWrite(fid, 'Growth info :');

% mu
dispAndWrite(fid, [' * From frame ' num2str(p.fitRange(1)) ' to ' num2str(p.fitRange(end)) ' gave mu :']);
dispAndWrite(fid, sprintf(['   - %1.2f as determined by number of cells'], muNrCells));
dispAndWrite(fid, sprintf(['   - %1.2f as determined by length of cells'], muLength));
dispAndWrite(fid, sprintf(['   - %1.2f as determined by volume of cells'], muVol));

% doublings
dispAndWrite(fid, [' * From frame ' num2str(p.fitRange(1)) ' to ' num2str(p.fitRange(end)) ' were doublings :']);
doubling = log2( NrCellsInFrame( p.fitRange(end) ) / NrCellsInFrame( p.fitRange(1) ) );
dispAndWrite(fid, sprintf(['   - %2.2f as determined by number of cells'], doubling));
doubling = log2( length_sum( p.fitRange(end) ) / length_sum( p.fitRange(1) ) );
dispAndWrite(fid, sprintf(['   - %2.2f as determined by length of cells'], doubling));
doubling = log2( area_sum( p.fitRange(end) ) / area_sum( p.fitRange(1) ) );
dispAndWrite(fid, sprintf(['   - %2.2f as determined by volume of cells'], doubling));
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Close file to write results to
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
fclose(fid);
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time_mins = convert_2_minutes(datetime);
s = sort([datetime(find(datetime>0))]);
firstmoment = s(1);
for i = 1:length(datetime),
  if (datetime(i)>0)
    tdiff = datetime(i) - firstmoment;
    time_mins(i) = DJK_getMinutesFromTimestamp(tdiff);
  else 
    time_mins(i)=0;
  end
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
