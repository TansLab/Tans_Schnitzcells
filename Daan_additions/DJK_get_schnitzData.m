% function schnitzData = DJK_get_schnitzData(p, s, timeField, varargin)
%
% DJK_get_schnitzData takes schnitzcells and returns selected data such
% that it can be plotted by DJK_plot_scatterColor. Mainly usefull as it can
% return multiple datapoints per schnitz (otherwise not plottable in
% DJK_plot_scatterColor). Also can select only datapoints within a
% particular timeFrame, so that only some points from a schnitz are used.
%
% A timefield needs to be provided, which will be added (together with
% schnitzNr) to the output schnitzData. This can be eg. av_time, Y_time &
% time.
%
% Data fields that are returned can be set using 'dataFields'. Data points
% are selected in two ways: whether a schnitz may be used or not, can be
% set with the field 'useForPlot' (only schnitzes that are 1 there, will be
% used). If 'useForPlot' is not set, all schnitzes with data will be used.
% If only data between certain timepoints must be used, this can be set
% using 'fitTime'. In this case it could be usefull to set timeField, in
% case less data values per schnitz than maximum.
%
% In general, the value at each dataField is returned for each timeField
% time point. If the average, first, last, min or max needs to returned,
% this can be accomplished by using 'useArrayData'. 'useArrayData' =
%   @(x,idx) mean(x)
%   @(x,idx) x(end)
%   @(x,idx) x(idx) (default)
%
% OUTPUT
% 'schnitzData'       cell structure with selected data
%
% REQUIRED ARGUMENTS:
% 'p'    
% 's'                 schnitzcells with useForPlot
% 'timeField'
%
% OPTIONAL ARGUMENTS:
% 'dataFields'        fields to be stored in schnitzData
%                     default: {'Y6_mean' 'muP13_fitNew'}
%
% 'fitTime'=[10 100]  time between datapoints in branches must lie
%                     (default: all timepoints)
%
% 'useArrayData'      
%

function schnitzData = DJK_get_schnitzData(p, s, timeField, varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 3; functionName = 'DJK_get_schnitzData';

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
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(s(1),char(timeField))
  disp(['timeField ' timeField ' does not exist. Exiting...!']);
  return;
end

% If not provided, use standard dataFields
if ~existfield(p, 'dataFields')
  p.dataFields = {'Y6_mean' 'muP13_fitNew'};
end

for i = 1:length(p.dataFields)
  field = char(p.dataFields(i));
  if ~existfield(s(1),field)
    disp(['Field ' field ' does not exist. Exiting...!']);
    return;
  end
end

if ~existfield(p,'fitTime')
  p.fitTime = [min([s.(timeField)]) max([s.(timeField)])];
end

if ~existfield(p,'useArrayData')
  p.useArrayData = @(x,idx) x(idx);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SOME PREPARATION FOR OUTPUT
%--------------------------------------------------------------------------
% schnitzData will contain data
schnitzData = struct;

% useAllcells=1 when useForPlot is not set
useAllcells = ~existfield(s(1),'useForPlot');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES AND GET DATA
%--------------------------------------------------------------------------
dataPointNr = 0;

% loop over schnitzes, from last to first
for schnitzNr = 1:length(s)

  % only consider schnitz if: 
  % * selected
  % * has data for the timeField
  if (useAllcells | s(schnitzNr).useForPlot)

    % loop over each timepoint for this schnitz
    for idx = 1:length(s(schnitzNr).(timeField))
      
      % check whether inside fitTime
      if s(schnitzNr).(timeField)(idx) >= p.fitTime(1) & s(schnitzNr).(timeField)(idx) <= p.fitTime(2) 

        dataPointNr = dataPointNr + 1;

        % add schnitzNrs
        schnitzData(dataPointNr).schnitzNr = schnitzNr;

        % add timeField
        schnitzData(dataPointNr).(timeField) = s(schnitzNr).(timeField)(idx);

        % add data fields
        for i = 1:length(p.dataFields)
          field = char(p.dataFields(i));
          if idx > length(s(schnitzNr).(field))
            schnitzData(dataPointNr).(field) = NaN;
          else
            schnitzData(dataPointNr).(field) = p.useArrayData( s(schnitzNr).(field), idx);
          end
        end

      end
    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% GET MEAN NOISE FIELDS
%--------------------------------------------------------------------------
timepoints = unique([schnitzData.(timeField)]);

% loop over timepoints
for t = 1:length(timepoints)

  idx = find([schnitzData.(timeField)] == timepoints(t)); % data points of this time point

  % loop over data fields
  for f = 1:length(p.dataFields)
    field = char(p.dataFields(f));
    noisefield = ['noise_' field];
    
    av = mean([schnitzData(idx).(field)]);

    % loop over data points
    for i = 1:length(idx)
      schnitzData(idx(i)).(noisefield) = schnitzData(idx(i)).(field) - av;
    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SHOW SOME OUTPUT
%--------------------------------------------------------------------------
disp(['---------------------------------------------------------------']);
disp(['Returned schnitzData struct with :']);
disp([' * ' num2str(length(schnitzData)) ' data points']);
disp([' * in between ' num2str(round(p.fitTime(1))) ' and ' num2str(round(p.fitTime(2))) ' mins.']);
disp([' * the following fields :']);
disp(['       - schnitzNr']);
disp(['       - ' timeField]);
for i = 1:length(p.dataFields)
  field = char(p.dataFields(i));
  disp(['       - ' field]);
end
disp(['---------------------------------------------------------------']);
%--------------------------------------------------------------------------

