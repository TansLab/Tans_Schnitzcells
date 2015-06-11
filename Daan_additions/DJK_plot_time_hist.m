% DJK_plot_time_hist
% 
%
% OUTPUT
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzData'  
% 'field'           	field from schnitzData used for histogram
% 'background'        autofluor
%
% OPTIONAL ARGUMENTS:
% 'onScreen' = 0      automatically save and close images
%              1      will ask to save and close images (default)
%              2      will not show or save images
% 'DJK_saveDir'       Directory where images will be saved. Default:
%                     "p.analysisDir 'schnitzcells\'"
% 'selectionName'
% 'binCenters'        default: [0:40]
% 'timeField'         default: 'Y_time'

function [mean_data noise_data data] = DJK_plot_time_hist(p, schnitzData, field, background, varargin);

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_plot_time_hist';

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
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(p,'timeField')
  p.timeField = 'Y_time';
end
if ~existfield(p,'binCenters')
  p.binCenters = [0:1:200];
end
% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 1;
end
% Just in case it has not been set yet
if ~existfield(p,'selectionName ')
  p.selectionName = '';
end
% Just in case it has not been set yet
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzcells' filesep];
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
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get data to plot
%--------------------------------------------------------------------------
timepoints = unique([schnitzData.(p.timeField)]);

% first seperate time points
for i = 1:length(timepoints)
  t = timepoints(i);
  idx = find([schnitzData.(p.timeField)] == t);
  data_hist(i,:) = hist([schnitzData(idx).(field)]-background, p.binCenters);
end

% now all data
binCenters_all_data = [min(p.binCenters): (max(p.binCenters) - min(p.binCenters) + 1) /  (4*length(p.binCenters)) : max(p.binCenters)];

all_data = [schnitzData.(field)]-background;
all_data_hist = hist(all_data, binCenters_all_data);
% normalize
all_data_hist = all_data_hist .* (max(max(data_hist)) / max(all_data_hist));

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------
% Make Figure Name
figureName = ['TimeHist ___ ' field];
figureFileName = ['TimeHist___' field];
xlabelName = field;

% make figure with full screen size
scrsz = get(0, 'ScreenSize');
fig1 = figure('Position', [126 scrsz(4)-450 scrsz(3)-125 scrsz(4)-450], 'Name', figureName, 'visible','off');
hold on;

% plot histograms
map = colorGray(length(timepoints));
for i = 1:length(timepoints)
  plot(p.binCenters, data_hist(i,:),'Color',map(i,:));
end

% label axes
xlabel(xlabelName,'interpreter','none','FontWeight','bold','FontSize',12); % interpreter to avoid problem with underscores
ylabel('#','interpreter','none','FontWeight','bold','FontSize',12);

% Add title
title([p.movieDate ' ' p.movieName ' for ' num2str(length(schnitzData)) ' schnitzData points in folder: schnitzcells\\' p.selectionName], ...
      'interpreter','none','FontWeight','bold','FontSize',12);

% Printout average, stdev and # of datapoints
text(0.8,0.975,['# points : ' num2str(length(all_data))],'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
text(0.8,0.925,['average  : ' num2str(mean(all_data))],'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
text(0.8,0.875,['stdev    : ' num2str(std(all_data))],'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
text(0.8,0.825,['noise    : ' num2str(round(100*std(all_data)/mean(all_data))) '%'],'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);

% plot hist all data
plot(binCenters_all_data, all_data_hist, '--', 'Color', 'k', 'LineWidth',2);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% ADD COLOR BAR
%--------------------------------------------------------------------------
p.colorNormalize = [min(timepoints) max(timepoints)];

colormap(map);
cBar = colorbar('EastOutside'); %, 'peer',axes_handle

% lowest value
labelText = {p.colorNormalize(1)}; 
labelHeight = [1];

% min color value
minColorHeight = round( DJK_scaleRange(min(timepoints), [p.colorNormalize(1) p.colorNormalize(2)], [1 1+length(timepoints)]) );
if minColorHeight>1 & minColorHeight<1+length(timepoints)
  labelText = [labelText 'min']; 
  labelHeight = [labelHeight minColorHeight]; 
end  

% max color value
maxColorHeight = round( DJK_scaleRange(max(timepoints), [p.colorNormalize(1) p.colorNormalize(2)], [1 1+length(timepoints)]) );
if maxColorHeight>1 & maxColorHeight>minColorHeight & maxColorHeight<1+length(timepoints)
  labelText = [labelText 'max']; 
  labelHeight = [labelHeight maxColorHeight]; 
end  

% highest value
labelText = [labelText p.colorNormalize(2)]; 
labelHeight = [labelHeight 1+length(timepoints)]; 

% set Ticks
set(cBar, 'YTickLabel', labelText, 'YTick', labelHeight);

% set Label at color bar
set(get(cBar,'YLabel'),'String',p.timeField,'interpreter','none','FontWeight','bold','FontSize',12);

% reduce width colorbar
graphPos = get(gca,'Position');
cBarPos = get(cBar,'Position');
cBarPos(3) = 0.5*cBarPos(3);
set(cBar,'Position', cBarPos);
set(gca,'Position', graphPos); % must do this, else wrong
%--------------------------------------------------------------------------

       
%--------------------------------------------------------------------------
% ASKING TO SAVE FIGURE
%--------------------------------------------------------------------------
% Ask to save the figure
if p.onScreen == 1
  set(fig1,'visible','on');
  saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
  pause(0.2);
elseif p.onScreen == 2
  saveFigInput='No';
else
  saveFigInput='Yes and Close';
end
if (upper(saveFigInput(1))=='Y')
    saveSameSize(fig1,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');
    if (strcmp(saveFigInput,'Yes and Close'))
        close(fig1);
        pause(0.2);
    end
  disp([' * Saved plot in ' figureFileName '.png']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PREPARE OUTPUT
%--------------------------------------------------------------------------
mean_data = mean(all_data);
noise_data = std(all_data)/mean(all_data);
data = all_data;
%--------------------------------------------------------------------------
