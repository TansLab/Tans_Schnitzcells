% DJK_plot_avColonyOverTime plots the average of particular data in
% schnitzcells over time, one data value per timepoint (frame). 
%
% If schnitzcells has the field 'useForPlot', only schnitzes that are 1
% there, will be used, otherwise all cells with data will be plotted. 
%
% If the field in schnitzcells to be plotted has only 1 datapoint per
% schnitz, this one is used for all of its timepoints.
%
% The average data value (and variance) is returned (with fitTime a
% particular part can be selected).
%
% Uses formula: mean     = sum(x) / #x
%               variance = sum( (x - mean)^2 )  / #x
%
% Old formula:  variance = (sum(x^2) - #x*mean) / (#x-1)
%
% For weighted mean / variance (using 'weightField'='lengthNew') uses formula:
%               mean     = sum(x * x_length/sum_x_length)
%               variance = sum( (x - mean)^2 )  / #
%
% OUTPUT
% 'fit_mean'        mean of data during fitTime TODO.......................
% 'fit_variance'    variance of data during fitTime
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzcells'    schnitzcells structure, where useForPlot must be set
% 'field'           field from schnitzcells used for Y axis
%
% OPTIONAL ARGUMENTS:
% 'manualRange'     Allows to analyze a subset of frames
% 'onScreen' = 0    automatically save and close image
%              1    will ask to save and close image (default)
% 'DJK_saveDir'     Directory where images will be saved.
%                   defaukl: "p.analysisDir 'schnitzcells\'"
% 'selectionName'
% 'xlim'
% 'ylim'
% 'timeField'       default: 'time', but can also be 'Y_time'
% 'fitTime'         used to return average of data. default: all data points
% 'weightField'     Returns weighted average & variance. default: 0 -> not weighted
% 'weightFrame'=1   For returning fit_mean use corrected version where each
%                   frame has same weight (default: 0)
%
% Hard-coded parameters
% FONTSIZE          Sets plotting fontsize

function [fit_mean, fit_variance] = DJK_plot_avColonyOverTime(p, schnitzcells, field, varargin);

FONTSIZE = 12;

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 3; functionName = 'DJK_plot_avColonyOverTime';

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
% If explicit manualRange is not given, take all in schnitzcells
if ~existfield(p,'manualRange')
  p.manualRange = unique( [schnitzcells.frame_nrs] ); % MW fix 2014/11
end
% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 1;
end
% 
if ~existfield(p,'timeField')
  p.timeField = 'time';
end
if ~existfield(p,'weightField')
  p.weightField = 'length_fitNew';
end
if ~existfield(p,'weightFrame')
  p.weightFrame = 0;
end
if ~existfield(p,'selectionName ')
  p.selectionName = '';
end
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
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Not plotting!');
  return;
end
if ~existfield(schnitzcells(1),field)
  disp(['Field ' field ' does not exist. Not plotting!']);
  return;
end
if p.weightField ~= 0
  if ~existfield(schnitzcells(1),p.weightField)
    disp(['weightField ' p.weightField ' does not exist. Not weighing!']);
    p.weightField = 0;
  else
    disp(['Weighting with field ' p.weightField]);
  end
end

% check whether field has only 1 value per schnitz
if length([schnitzcells.(field)]) <= length([schnitzcells.P])
  onlyOneDataPoint = 1;
  disp(['Field ' field ' has only 1 datapoint per schnitz!']);
else
  onlyOneDataPoint = 0;
end

% check whether field has less value per schnitz than timepoints
if length([schnitzcells.(field)]) < length([schnitzcells.(p.timeField)])
  disp(['Field ' field ' has less datapoints per schnitz than timepoints!']);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parse extra input arguments
%--------------------------------------------------------------------------
if ~existfield(p,'fitTime')
  allTime = unique([schnitzcells.(p.timeField)]);
  p.fitTime = [allTime(1) allTime(end)];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Collect all relevant data
%--------------------------------------------------------------------------
data = struct;
nr_frames = length(p.manualRange);
for i = 1:nr_frames
  data(i).frame = p.manualRange(i);
  data(i).time = -1;
  data(i).x = [];
  data(i).weight = [];
end

useAllcells = ~existfield(schnitzcells(1),'useForPlot');

% loop over cells
for cell = 1:length(schnitzcells)
  % check whether cell will be used for plotting
  if useAllcells | schnitzcells(cell).useForPlot
    % check whether cell has relevant data
    if length(schnitzcells(cell).(field)) > 0
      % loop over frames where cell lives
      for age = 1:length(schnitzcells(cell).frame_nrs)
          
        fr = schnitzcells(cell).frame_nrs(age); % fix mw 2014/11
        i = find( p.manualRange == fr);

        if i
          data(i).time = schnitzcells(cell).(p.timeField)(age);

          if onlyOneDataPoint
            x = [data(i).x schnitzcells(cell).(field)(1)];
          else
            x = [data(i).x schnitzcells(cell).(field)(age)];
          end
          
          if ~isnan(x)
            data(i).x = x;

            if p.weightField ~= 0
              if onlyOneDataPoint
                data(i).weight = [data(i).weight schnitzcells(cell).(p.weightField)(1)];
              else
                data(i).weight = [data(i).weight schnitzcells(cell).(p.weightField)(age)];
              end
            else
              data(i).weight = [data(i).weight 1];
            end
          end
          
        end
      end
    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get mean & variance of each frame data to plot
%--------------------------------------------------------------------------
for i = 1:nr_frames
  if length(data(i).x > 0)
    data(i).mean = sum(data(i).x .* data(i).weight) / sum(data(i).weight);
    data(i).median = median(data(i).x); %MW
    data(i).variance = sum( (data(i).x - data(i).mean).^2 .* data(i).weight) / sum(data(i).weight);
  else
    data(i).mean = NaN;
    data(i).median = NaN; % MW
    data(i).variance = NaN;
  end
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get data to plot
%--------------------------------------------------------------------------
data_frames         = [];
data_time           = [];
data_field_mean     = [];
data_field_median   = []; % MW
data_field_variance = [];
for i = 1:nr_frames
  if ~isnan(data(i).mean) & length(data(i).x)>0
    data_frames         = [data_frames data(i).frame];
    data_time           = [data_time data(i).time];
    data_field_mean     = [data_field_mean data(i).mean];
    data_field_median   = [data_field_median data(i).median]; % MW
    data_field_variance = [data_field_variance data(i).variance];
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Determine fit value
%--------------------------------------------------------------------------
fit_x       = [];
fit_weight  = [];
fit_weight_frame_corrected  = [];

for i = 1:nr_frames
  if data(i).time >= p.fitTime(1) & data(i).time <= p.fitTime(2)
    fit_x       = [fit_x data(i).x];
    fit_weight  = [fit_weight data(i).weight];
    fit_weight_frame_corrected  = [fit_weight_frame_corrected (data(i).weight/sum(data(i).weight))];
  end
end

% determine mean & variance letting each cell contribute equally
fit_mean = sum(fit_x .* fit_weight) / sum(fit_weight);
fit_variance = sum( (fit_x - fit_mean).^2 .* fit_weight) / sum(fit_weight);

% determine mean & variance letting each frame contribute equally
fit_mean_frame_corrected = sum(fit_x .* fit_weight_frame_corrected) / sum(fit_weight_frame_corrected);
fit_variance_frame_corrected = sum( (fit_x - fit_mean_frame_corrected).^2 .* fit_weight_frame_corrected) / sum(fit_weight_frame_corrected);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------
% Make Figure Name
figureName = ['av_' field ' ___time']; % DJK 090319 was field ' ___ ' field
figureFileName = [p.movieName '_av_' field '___time' ]; % DJK 090319 was field ' ___ ' field
figureTitle = [p.movieDate ' ' p.movieName ' for ' num2str(length([data.x])) ' schnitzes in folder: schnitzcells\\' p.selectionName];
figureXlabel = ['time (mins)'];
figureYlabel = ['average ' field];

% Adjust if weighted
if p.weightField ~= 0
  figureName = [figureName '___weighted_' p.weightField]; 
  figureFileName = [figureFileName '___weighted_' p.weightField];
  figureYlabel = [figureYlabel ' weighted with ' p.weightField];
end

% make figure with full screen size
%{
scrsz = get(0, 'ScreenSize');
pos1=[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.8 scrsz(4)*0.7];
fig1 = figure('Position', pos1, 'Name', figureName, 'visible','off');
%}
fig1=figure();

hold on;

% plot standard deviation as filled area
std_data_time = [data_time data_time(end:-1:1)]; % DJK 090319 was data_frames
std_data_std = [data_field_mean-sqrt(data_field_variance) data_field_mean(end:-1:1)+sqrt(data_field_variance(end:-1:1))];
fill(std_data_time,std_data_std,'r','EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8]);

% plot colony-median as function frames (i.e. time) - MW
line( 'Xdata',data_time, ... % DJK 090319 was data_frames
      'Ydata',data_field_median,... 
      'LineStyle','none', ...
      'Color','b', ...
      'Marker','o' );

% plot colony-mean as function frames (i.e. time)
line( 'Xdata',data_time, ... % DJK 090319 was data_frames
      'Ydata',data_field_mean,... 
      'LineStyle','none', ...
      'Color','k', ...
      'Marker','s' );

% Plot overall mean line
if fit_x
 hold on;
  %line( 'Xdata',p.fitTime, ...
  line( 'Xdata',[min(data_time), max(data_time)], ...
        'Ydata',[fit_mean fit_mean], ...
        'LineStyle',':', ...
        'LineWidth',4, ...
        'Color','r', ...
        'Marker','none');
  line1= ['Fit from (min): ' sprintf('%0.0f', p.fitTime(1)) ' to ' sprintf('%0.0f',p.fitTime(2))];
  line2= ['Fitted mean & variance: ' sprintf('%0.2f',fit_mean) ' / ' sprintf('%0.2f',fit_variance)];
  line3= ['Fitted mean & variance: ' sprintf('%0.2f',fit_mean_frame_corrected) ' / ' sprintf('%0.2f',fit_variance_frame_corrected)];
  line2a= [' (cells equally contributing)'];
  line3a= [' (frames equally contributing)'];
  figureTitle = [figureTitle,10,line1,10,line2,line2a,10,line3,line3a];
    
end

% label axes
xlabel(figureXlabel,'interpreter','none','FontWeight','bold','FontSize',FONTSIZE); % interpreter to avoid problem with underscores
ylabel(figureYlabel,'interpreter','none','FontWeight','bold','FontSize',FONTSIZE);
  
% Add title
title(figureTitle,'interpreter','none','FontWeight','bold','FontSize',FONTSIZE);
      %'Units', 'normalized', 'Position', [1 1], 'HorizontalAlignment', 'right');


set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal'); set(gca,'FontSize',FONTSIZE)



% in case xlim or ylim
if existfield(p,'xlim '), xlim(p.xlim); end
if existfield(p,'ylim '), ylim(p.ylim); end


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


%--------------------------------------------------------------------------
% DATA TO BE RETURNED
%--------------------------------------------------------------------------
if p.weightFrame
  fit_mean = fit_mean_frame_corrected;
  fit_variance = fit_variance_frame_corrected;
end
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function number_string = DJK_setDecimalPlaces(number, decimal_places);
number_string = sprintf(['%1.' num2str(decimal_places) 'f'], number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















% %--------------------------------------------------------------------------
% % Get mean & variance of each frame data to plot
% %--------------------------------------------------------------------------
% useAllcells = ~existfield(schnitzcells(1),'useForPlot');
% all_frames = unique( [schnitzcells.frames] -1);
% 
% % initialize alldata
% for fr = all_frames;
%   alldata_field_sum(fr) = 0;
%   alldata_cellsContributing(fr) = 0;
%   alldata_field_sum2(fr) = 0;
% 
%   alldata_field_weighted_sum(fr) = 0;
%   alldata_field_weight(fr) = 0;
% end
% 
% nrCells = 0;
% % loop over cells to calc sum and # in each frame
% for cell = 1:length(schnitzcells)
% 
%   % check whether cell will be used for plotting
%   if useAllcells | schnitzcells(cell).useForPlot
%     nrCells = nrCells + 1;
% 
%     % check whether cell has relevant data
%     if length(schnitzcells(cell).(field)) > 0
% 
%       % loop over frames where cell lives
%       for age = 1:length(schnitzcells(cell).frames)
%         
%         fr = schnitzcells(cell).frames(age)-1;
%         alldata_cellsContributing(fr) = alldata_cellsContributing(fr) + 1;
%         alldata_time(fr) = schnitzcells(cell).(p.timeField)(age);
%         
%         if onlyOneDataPoint
%           alldata_field_sum(fr) = alldata_field_sum(fr) + schnitzcells(cell).(field)(1);
%           
%           if p.weightField ~= 0
%             alldata_field_weighted_sum(fr) = alldata_field_weighted_sum(fr) + schnitzcells(cell).(field)(1)*schnitzcells(cell).(p.weightField)(1);
%             alldata_field_weight(fr) = alldata_field_weight(fr) + schnitzcells(cell).(p.weightField)(1);
%           end
%         else
%           alldata_field_sum(fr) = alldata_field_sum(fr) + schnitzcells(cell).(field)(age);
%           
%           if p.weightField ~= 0
%             alldata_field_weighted_sum(fr) = alldata_field_weighted_sum(fr) + schnitzcells(cell).(field)(age)*schnitzcells(cell).(p.weightField)(age);
%             alldata_field_weight(fr) = alldata_field_weight(fr) + schnitzcells(cell).(p.weightField)(age);
%           end
%         end
%       end
%     end
%   end
% end
% 
% % calc mean for each frame
% for fr = all_frames;
%   alldata_mean(fr) = alldata_field_sum(fr) / alldata_cellsContributing(fr);
%   
%   if p.weightField ~= 0
%     alldata_weighted_mean(fr) = alldata_field_weighted_sum(fr) / alldata_field_weight(fr);
%   end
% end
% 
% % loop over cells to calc variance in each frame
% for cell = 1:length(schnitzcells)
%   if useAllcells | schnitzcells(cell).useForPlot
%     if length(schnitzcells(cell).(field)) > 0
%       for age = 1:length(schnitzcells(cell).frames)
%         fr = schnitzcells(cell).frames(age)-1;
%         if onlyOneDataPoint
%           alldata_field_sum2(fr) = alldata_field_sum2(fr) + (schnitzcells(cell).(field)(1) - alldata_mean(fr))^2;
%         else
%           alldata_field_sum2(fr) = alldata_field_sum2(fr) + (schnitzcells(cell).(field)(age) - alldata_mean(fr))^2;
%         end
%       end
%     end
%   end
% end
% 
% % calc variance for each frame
% for fr = all_frames;
%   if alldata_cellsContributing(fr) < 1
%     alldata_variance(fr) = nan;
%   else
%     alldata_variance(fr) = alldata_field_sum2(fr) / alldata_cellsContributing(fr);
%   end
% end
% %--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % Determine fit value
% %--------------------------------------------------------------------------
% % determine mean & variance letting each cell contribute equally
% fit_n     = 0;
% fit_sum   = 0;
% fit_sum2  = 0;
% for cell = 1:length(schnitzcells)
%   if useAllcells | schnitzcells(cell).useForPlot
%     if length(schnitzcells(cell).(field)) > 0
%       for age = 1:length(schnitzcells(cell).frames)
%         if schnitzcells(cell).time(age) >= p.fitTime(1) & schnitzcells(cell).time(age) <= p.fitTime(2)
%           fr = schnitzcells(cell).frames(age)-1;
%           fit_n = fit_n + 1;
%           if onlyOneDataPoint
%             fit_sum = fit_sum + schnitzcells(cell).(field)(1);
%           else
%             fit_sum = fit_sum + schnitzcells(cell).(field)(age);
%           end
%         end
%       end
%     end
%   end
% end
% fit_mean = fit_sum / fit_n;
% 
% if fit_n > 0
%   for cell = 1:length(schnitzcells)
%     if useAllcells | schnitzcells(cell).useForPlot
%       if length(schnitzcells(cell).(field)) > 0
%         for age = 1:length(schnitzcells(cell).frames)
%           if schnitzcells(cell).time(age) >= p.fitTime(1) & schnitzcells(cell).time(age) <= p.fitTime(2)
%             fr = schnitzcells(cell).frames(age)-1;
%             if onlyOneDataPoint
%               fit_sum2 = fit_sum2 + (schnitzcells(cell).(field)(1) - fit_mean)^2;
%             else
%               fit_sum2 = fit_sum2 + (schnitzcells(cell).(field)(age) - fit_mean)^2;
%             end
%           end
%         end
%       end
%     end
%   end
%   fit_variance = fit_sum2 / fit_n;
% else
%   fit_variance = NaN;
% end
% 
% % determine mean & variance letting each frame contribute equally
% % 2DO
% %------------------------------------------------------------------------
% --


% %--------------------------------------------------------------------------
% % Get data to plot
% %--------------------------------------------------------------------------
% % determine which values in alldata are in manualRange & have data for field
% count = 0;
% idx = [];
% for fr = p.manualRange
%   i = find( all_frames == fr);
%   if alldata_cellsContributing(fr) > 0 & ~isnan(alldata_field_sum(fr))
%     count = count+1;
%     idx(count) = i;
%   end
% end
% 
% % get data to plot
% data_frames         = p.manualRange;
% data_time           = alldata_time( p.manualRange(idx) );
% data_field_mean     = alldata_mean( p.manualRange(idx) );
% data_field_variance = alldata_variance( p.manualRange(idx) );
% %------------------------------------------------------------------------
% --
