% DJK_xcorr_branches can perform different kind of
% correlations: cross / auto, normalized or not, corrected or not
%
% uses Formula:   r * interval = tau
%                 Sxy(r) = n=1_SUM_N-r ( X(t) * Y(t+r) ) / (N-|r|)
%
% normalize:      Rxy(r) = Sxy(r) / ( sqrt(Sxx(0)) * sqrt(Syy(0)) )
%                 when mean(x)=mean(y)=0 this is the same as stdev_x * stdev_y
%
% not corrected:  in Sxy(r) don't divide by (N-|r|) -> Nature439_608 way
%                 for larger delay times the autocorrelation will always decrease
%
% autocorr:       provide dataField twice (as field1 & field2)
%
% REQUIRED ARGUMENTS:
% 'p'
% 'branchData'
% 'fieldX'        field X of which crosscorrelation should be calculated
% 'fieldY'        field Y of which crosscorrelation should be calculated
%
% OPTIONAL ARGUMENTS:
% 'timeField'     field used as time and stored in branchData
%                 default: 'Y_time'
% 'normalize'=0   no normalization (default: 1)
% 'correct'=0     no correction for less data (default: 1)
% 'onScreen' = 0  will automatically save and close images, default: ask.
% 'DJK_saveDir'   directory where images will be saved. 
%                 default: [p.analysisDir 'branches\']
% 'selectionName' added as subfolder to DJK_saveDir
% 'xlim'  
% 'ylim'
% 'samplingError' how many minutes can time be different from constant interval
%                 default: 0.3
% 'weighing'=0    Does not use the weighing feature (default=1)

function branchData = DJK_xcorr_branches(p, branchData, fieldX, fieldY, varargin);

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_xcorr_branches';

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


% --------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
% --------------------------------------------------------------------------
if ~existfield(p,'timeField')
  p.timeField = 'Y_time';
end

% check whether fields exist
if ~existfield(branchData(1),fieldX)
  disp(['Field ' fieldX ' does not exist. Exiting...!']); return;
end
if ~existfield(branchData(1),fieldY)
  disp(['Field ' fieldY ' does not exist. Exiting...!']); return;
end
if ~existfield(branchData(1), p.timeField)
  disp(['Field ' p.timeField ' does not exist. Exiting...!']); return;
end

% if too little data, exit
if (length(branchData(1).(fieldX)) < 2), error('Not enough data'); end

% settings
if ~existfield(p,'normalize')
  p.normalize = 1;
end
if ~existfield(p,'correct')
  p.correct = 1;
end
if ~existfield(p,'weighing')
  p.weighing = 1;
end

if ~existfield(p,'ylim') & p.normalize
  p.ylim = [-1 1];
end
if ~existfield(p,'samplingError')
  p.samplingError = 0.3;
end

% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 1;
end

% Save directory
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'branches' filesep];
end
if ~existfield(p,'selectionName ')
  p.selectionName = '';
end
if length(p.selectionName) > 0
  p.DJK_saveDir = [p.DJK_saveDir p.selectionName filesep];
end  
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Check sampling interval
% --------------------------------------------------------------------------
% check whether constant sampling in timefield
timeData = branchData(1).(p.timeField);
timeData = timeData-timeData(1);
fit = polyfit([1:length(timeData)],timeData,1);
interval = fit(1);
fittedTimeData = fit(2) + fit(1)*[1:length(timeData)];
interval_error = 0;
for i = 1:length(timeData)
  if abs(timeData(i) - fittedTimeData(i)) > p.samplingError
    disp(['Watch out! timepoints might not be evenly spaced! for i=' num2str(i) ' : timeData=' num2str(timeData(i)) ' & fittedTimeData=' num2str(fittedTimeData(i))]);
    interval_error = 1;
  end
end
if interval_error & p.onScreen
  figure; plot(timeData,'bo'); hold on; plot(fittedTimeData,'k-');
end
disp(['Time interval is ' num2str(interval) ' (mins)']);
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Prepare some stuff
% --------------------------------------------------------------------------
% crosscorrelation for each independent branch will be saved in branchData in this field
p.targetField = ['xCorr_' fieldX '_' fieldY];

% check whether branches have same length
max_br_length = 0;
min_br_length  = realmax('double');
for br = 1:length(branchData)
  max_br_length = max(max_br_length, length(branchData(br).(p.timeField)));
  min_br_length = min(min_br_length, length(branchData(br).(p.timeField)));
end
if (max_br_length == min_br_length)
  disp(['All branches have same length: ' num2str(min_br_length) ' timepoints.']);
else
  disp(['Branches have different lengths: from ' num2str(min_br_length) ' to ' num2str(max_br_length) ' timepoints.']);
end

% time for crosscorrelation
for br = 1:length(branchData)
  branchData(br).corr_time = interval * [-max_br_length+1:max_br_length-1];
end
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Calc of correlation for each branch
% --------------------------------------------------------------------------
for br = 1:length(branchData)
  X = branchData(br).(fieldX);
  Y = branchData(br).(fieldY);
  W = ones(size(X));

  % Get Sxy for this branch
  branchData(br).(p.targetField) = func_Sxy( X, Y, W, max_br_length, p.correct);
  
  % in case normalization: divide by sqrt(Sxx(0)) * sqrt(Syy(0))
  if p.normalize
    if p.correct
      branchData(br).(p.targetField) = branchData(br).(p.targetField) ./ ...
                                     sqrt( func_Sxyr_correct(X,X,W,0) * func_Sxyr_correct(Y,Y,W,0));
    else
      branchData(br).(p.targetField) = branchData(br).(p.targetField) ./ ...
                                     sqrt( func_Sxyr(X,X,W,0) * func_Sxyr(Y,Y,W,0));
    end
  end
end
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Calc of correlation of composite
% --------------------------------------------------------------------------
X = struct;
Y = struct;
W = struct;
for br = 1:length(branchData)
  X(br).data = branchData(br).(fieldX);
  Y(br).data = branchData(br).(fieldY);
  if p.weighing,  W(br).data = branchData(br).weight;
  else,           W(br).data = ones(size(X(br).data));
  end
end

% Get Sxy_composite for all branches
Sxy_composite = func_Sxy_composite( X, Y, W, max_br_length, p.correct);

% in case normalization: divide by sqrt(Sxx_composite(0)) * sqrt(Syy_composite(0))
if p.normalize
  if p.correct
    Sxy_composite = Sxy_composite ./ sqrt( func_Sxyr_composite_correct(X,X,W,0) * func_Sxyr_composite_correct(Y,Y,W,0));
  else
    Sxy_composite = Sxy_composite ./ sqrt( func_Sxyr_composite(X,X,W,0) * func_Sxyr_composite(Y,Y,W,0));
  end
end

% ELOWITZ - Gives same result as func_Sxyr_composite, for weighing=0 & correct=1 and sorted Branches
% K = struct;
% for br = 1:length(branchData), K(br).data = branchData(br).branchPoint; end
% Sxy_composite = func_Sxy_composite_correct_Elowitz( X, Y, W, K, max_br_length);
% 
% % ELOWITZ - in case normalization: divide by sqrt(Sxx_composite(0)) * sqrt(Syy_composite(0))
% if p.normalize
%   % always correct
%   Sxy_composite = Sxy_composite ./ sqrt( func_Sxyr_composite_correct_Elowitz(X,X,W,K,0) * func_Sxyr_composite_correct_Elowitz(Y,Y,W,K,0));
% end
% -------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------
% Make Figure Name
figureName = ['xCorr ___ ' fieldX '_' fieldY];
figureFileName = ['xCorr_' fieldX '_' fieldY];

if p.correct
  figureName = [figureName '_correct1'];
  figureFileName = [figureFileName '_correct1'];
end

scrsz = get(0, 'ScreenSize');
% fig1 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
fig1 = figure('Position', [126 scrsz(4)-150 scrsz(3)-125 scrsz(4)-150], 'Name', figureName, 'visible','off');

% use this colormap
map = DJK_hsv(length(branchData));

% loop over branches
for br = 1:length(branchData)
  plot(branchData(br).corr_time, branchData(br).(p.targetField), 'Color',map(br,:), 'LineWidth',1); %
  hold on;
end

% plot composite
plot(branchData(1).corr_time, Sxy_composite, 'Color', 'k', 'LineWidth',4); %
hold on;

ylabel(p.targetField,'interpreter','none','FontWeight','bold','FontSize',12);
xlabel(['time (mins)'],'interpreter','none','FontWeight','bold','FontSize',12);

% in case xlim or ylim
if existfield(p,'xlim ')
  xlim(p.xlim);
end
if existfield(p,'ylim ')
  ylim(p.ylim);
end

% plot vertical line at x=0
plot([0 0], [-0.9 0.9], 'k-', 'LineWidth',1);
hold on;

% plot horizontal line at y=0
plot(xlim, [0 0], 'k-', 'LineWidth',1);
hold on;

% Add title
title([p.movieDate ' ' p.movieName ' for ' num2str(length(branchData)) ' branches in folder: branches\\' p.selectionName],'interpreter','none','FontWeight','bold','FontSize',12);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% ASKING TO SAVE FIGURE
%--------------------------------------------------------------------------
% Ask to save the figure
if p.onScreen
  set(fig1,'visible','on');
  saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
  pause(0.2);
else
  saveFigInput='Yes and Close';
end
if (upper(saveFigInput(1))=='Y')
%     saveas(fig1,[p.DJK_saveDir figureFileName '.fig']);
    saveSameSize(fig1,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');
    if (strcmp(saveFigInput,'Yes and Close'))
        close(fig1);
        pause(0.2);
    end
  disp([' * Saved plot in ' figureFileName '.png']);
end
%--------------------------------------------------------------------------

% both weighing, do:
%   replace W(1+r:N) by (0.5*W(1+r:N) + 0.5*W(1:N-r))
%   replace W(br).data(1+r:N) by (0.5*W(br).data(1+r:N) + 0.5*W(br).data(1:N-r))
%   replace W(N:-1:-r+1) by (0.5*W(N:-1:-r+1) + 0.5*W(N+r:-1:1))
%   replace W(br).data(N:-1:-r+1) by (0.5*W(br).data(N:-1:-r+1) + 0.5*W(br).data(N+r:-1:1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I : Sxy(tau)
function Sxyr = func_Sxyr( X, ...   % data X
                           Y, ...   % data Y
                           W, ...   % weight
                           r)       % r = tau / interval
N = length(X); % length of this particular branch

% calculate Sxyr in case r is positive
if r >= 0
  Sxyr = sum( X(1:N-r) .* Y(1+r:N) .* W(1+r:N) );

% calculate Sxyr in case r is negative
else
  Sxyr = func_Sxyr( Y, X, W, -r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II : Sxy_correct(tau)
function Sxyr = func_Sxyr_correct( X, ...   % data X
                                   Y, ...   % data Y
                                   W, ...   % weight
                                   r)       % r = tau / interval
N = length(X); % length of this particular branch
Sxyr = func_Sxyr( X, Y, W, r) / sum( W(1+abs(r):N) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III : Sxy_composite(tau)
function Sxyr_composite = func_Sxyr_composite( X, ...   % struct with X in X.data
                                               Y, ...   % struct with Y in Y.data
                                               W, ...   % struct with weight in W.data
                                               r)       % r = tau / interval
Sxyr_composite = 0;
for br = 1:length(X)
  Sxyr_composite = Sxyr_composite + func_Sxyr(X(br).data,Y(br).data,W(br).data,r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV : Sxy_composite_correct(tau)
function Sxyr_composite = func_Sxyr_composite_correct( X, ...   % struct with X in X.data
                                                       Y, ...   % struct with Y in Y.data
                                                       W, ...   % struct with weight in W.data
                                                       r)       % r = tau / interval
sum_weight = 0;
for br = 1:length(X)
  N = length(W(br).data); % length of this particular branch
  sum_weight = sum_weight + sum( W(br).data(1+abs(r):N) );
end
Sxyr_composite = func_Sxyr_composite( X, Y, W, r) / sum_weight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V : Sxy_composite_correct_Elowitz(tau)
% Should be used with all Weights being 1, and all branches same length
function Sxyr_composite = func_Sxyr_composite_correct_Elowitz( X, ...   % struct with X in X.data
                                                               Y, ...   % struct with Y in Y.data
                                                               W, ...   % struct with weight in W.data
                                                               K, ...   % struct with branching point in K.data
                                                               r)       % r = tau / interval
sum_weight = 0;
subtract_Sxyr_composite = 0;
subtract_sum_weight = 0;
for br = 1:length(X)
  N = length(W(br).data); % length of this particular branch
  sum_weight = sum_weight + sum( W(br).data(1+abs(r):N) );

  branchPoint = K(br).data;
  subtract_Sxyr_composite = subtract_Sxyr_composite + func_Sxyr( X(br).data(1:branchPoint), Y(br).data(1:branchPoint), W(br).data(1:branchPoint), r);
  subtract_sum_weight = subtract_sum_weight + sum( W(br).data(1+abs(r):branchPoint) );
end
Sxyr_composite = (func_Sxyr_composite( X, Y, W, r) - subtract_Sxyr_composite) / (sum_weight - subtract_sum_weight);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sxy(all)
function Sxy = func_Sxy( X, ...              % data X
                         Y, ...              % data Y
                         W, ...              % data weight
                         max_br_length, ...  % length of longest branch
                         correct)            % correct for less datapoints?
% loop over r's
for r = -max_br_length+1:max_br_length-1
  if correct
    Sxy(max_br_length+r) = func_Sxyr_correct( X, Y, W, r);
  else
    Sxy(max_br_length+r) = func_Sxyr( X, Y, W, r);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sxy_composite(all)
function Sxy_composite = func_Sxy_composite( X, ...              % struct with X in X.data
                                             Y, ...              % struct with Y in Y.data
                                             W, ...              % struct with weight in W.data
                                             max_br_length, ...  % length of longest branch
                                             correct)            % correct for less datapoints?
% loop over r's
for r = -max_br_length+1:max_br_length-1
  if correct
    Sxy_composite(max_br_length+r) = func_Sxyr_composite_correct( X, Y, W, r);
  else
    Sxy_composite(max_br_length+r) = func_Sxyr_composite( X, Y, W, r);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sxy_composite_correct_Elowitz(all) % always correct
function Sxy_composite = func_Sxy_composite_correct_Elowitz( X, ...              % struct with X in X.data
                                                             Y, ...              % struct with Y in Y.data
                                                             W, ...              % struct with weight in W.data
                                                             K, ...              % struct with branchPoint in K.data
                                                             max_br_length), ...  % length of longest branch
% loop over r's
for r = -max_br_length+1:max_br_length-1
  % always correct
  Sxy_composite(max_br_length+r) = func_Sxyr_composite_correct_Elowitz( X, Y, W, K, r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
