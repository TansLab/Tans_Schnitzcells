function branchData = DJK_get_branches(p, s, varargin)
% function branchData = DJK_get_branches(p, s, varargin)
%
% DJK_get_branches takes schnitzcells and returns branch data of
% particular branches & defined data. 
%
% Data fields that are returned can be set using 'dataFields'. schnitzNrs
% datafield will always be added to branches. Of each field except the
% timefield (dataFields(1)) also a noise and norm version are returned.
% noise: mean of returned branches for this timepoint is subtracted. norm:
% mean of returned branches for all timepoints is subtracted.  
%
% Branches are selected in two ways: whether a schnitz may be used or not,
% can be set with the field 'useForPlot' (only schnitzes that are 1 there,
% will be used). If 'useForPlot' is not set, all schnitzes with data will be
% used. If only data between certain timepoints must be used, this can be
% set using 'fitTime'
%
% Note that only branched with maximum length are returned. Using
% 'sameLength'=0 can be specified that branches of unequal length are returned. 
%
% Some data might be used more than once (in different branches). This is
% counted in the weight field.
%
% 2DO: ranking
%
% OUTPUT
% 'branchData'        cell structure with branches
%
% REQUIRED ARGUMENTS:
% 'p'                 parameter struct, with 
%                       p.dataFields (default {'Y_time' 'Y_mean' 'mu_fitCoef3'};)
%                       p.singleSchnitz (optional)
%                       p.fitTime (optional)
%                       p.sameLength (optional, default 1): [not sure but I think if =1 takes into account only branches that run from fitTime(1) to fitTime(2), i.e. are "complete"]
%                       (also stores parameters given by varargin)
%
% 's'                 schnitzcells with useForPlot
%
% OPTIONAL ARGUMENTS:
% 'dataFields'        fields to be stored in branchData
%                     default: {'Y_time' 'Y_mean' 'mu_fitCoef3'}
%                     note: should have same # values!!!! Put time in first one!
%                     noiseField will be added & renamed: 'Y_mean' becomes 'noise_Y_mean'
%                     normField will be added & renamed: 'Y_mean' becomes 'norm_Y_mean'
%
% 'sameLength'=0      does not make branches the same length
%
% 'singleSchnitz'=1   makes branch length only 1 schnitz long & not same length (default:0)
%
% 'fitTime'=[10 100]  time between datapoints in branches must lie
%                     (default: all timepoints)
%


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'DJK_get_branches';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
  warning(errorMessage);
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
% If not provided, use standard dataFields
if ~existfield(p, 'dataFields')
  p.dataFields = {'Y_time' 'Y_mean' 'mu_fitCoef3'};
end

for i = 1:length(p.dataFields)
  field = char(p.dataFields(i));
  if ~existfield(s(1),field)
    disp(['Field ' field ' does not exist. Exiting...!']);
    return;
  end
end

% first dataField should contain time (timeField)
timeField = char(p.dataFields(1));

if ~existfield(p,'singleSchnitz')
  p.singleSchnitz = 0;
end

if ~existfield(p,'sameLength')
  if p.singleSchnitz, p.sameLength = 0;
  else,               p.sameLength = 1; end
end

if ~existfield(p,'fitTime')
  p.fitTime = [min([s.(timeField)]) max([s.(timeField)])];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SOME PREPARATION FOR OUTPUT
%--------------------------------------------------------------------------
% remembers whether data for this schnitz has been used already
for i = 1:length(s)
  s(i).schnitzUsed = 0;
end

% branchData will contain data
branchData = struct;

% useAllcells=1 when useForPlot is not set
useAllcells = ~existfield(s(1),'useForPlot');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES AND GET DATA
%--------------------------------------------------------------------------
branchNr = 0;

% loop over schnitzes, from last to first
for schnitzNr = length(s):-1:1

  % only consider schnitz if: 
  % * selected 
  % * not incorporated in another branch
  if (useAllcells | s(schnitzNr).useForPlot) & ...
     (~s(schnitzNr).schnitzUsed) & ...
     (length(s(schnitzNr).(timeField)) > 0)

    % * first time point not outside p.fitTime (else branches can be double)
    if (s(schnitzNr).(timeField)(1) <= p.fitTime(2))
  
      % This schnitz will function as end schnitz of branch 'branchNr'
      branchNr = branchNr+1;
      done = 0;

      cur_schnitzNr = schnitzNr;
      data_temp = cell(1+length(p.dataFields),1);

      % now get data for this branch from end to beginning
      while ~done
        s(cur_schnitzNr).schnitzUsed = 1;
        temp = [];

        % add data fields
        for i = 1:length(p.dataFields)
          field = char(p.dataFields(i));
          temp = getfield(s(cur_schnitzNr),field);
          temp = temp(end:-1:1);
          data_temp{i+1,:} = cat(2,data_temp{i+1,:}, temp);
        end

        % add schnitzNrs
        data_temp{1,:} = cat(2, data_temp{1,:}, cur_schnitzNr*ones(1,length(temp)));

        % move to parent
        cur_schnitzNr = s(cur_schnitzNr).P;

        % check whether we are done
        done = (cur_schnitzNr <=0);
        if ~done
          if ~useAllcells & ~s(cur_schnitzNr).useForPlot
            done = 1;
          end
        end

        % in case of singleSchnitz also done after 1 loop
        if p.singleSchnitz, done = 1; end
      end

      % invert data_temp 
      for i = 1:length(data_temp)
        data_invers = data_temp{i,:};
        data_temp{i,:} = data_invers(end:-1:1);
      end

      % add to branchData
      branchData(branchNr).schnitzNrs = data_temp{1,:};
      for i = 1:length(p.dataFields)
        field = char(p.dataFields(i));
        branchData(branchNr).(field) = data_temp{i+1,:};
      end

    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER BRANCHES AND REMOVE DATA OUTSIDE TIME FRAME
%--------------------------------------------------------------------------
% DETERMINE TIME FRAME FOR ALL BRANCHES
all_time_data = sort([branchData.(timeField)]);
idx = find(all_time_data >= p.fitTime(1));
timeFrame(1) = all_time_data(idx(1));
idx = find(all_time_data <= p.fitTime(2));
timeFrame(2) = all_time_data(idx(end));

disp(['---------------------------------------------------------------']);
disp(['fitTime            : ' num2str(round(p.fitTime(1))) ' mins till ' num2str(round(p.fitTime(2))) ' mins']);
disp(['useForPlot time    : ' num2str(round(min(all_time_data))) ' mins till ' num2str(round(max(all_time_data))) ' mins']);
disp(['timeFrame          : ' num2str(round(timeFrame(1))) ' mins till ' num2str(round(timeFrame(2))) ' mins']);
disp(['---------------------------------------------------------------']);



% loop over branches
for branchNr = 1:length(branchData)
  min_idx = find(branchData(branchNr).(timeField)>=timeFrame(1));
  max_idx = find(branchData(branchNr).(timeField)<=timeFrame(2));
  if length(min_idx)==0 | length(max_idx)==0, disp('error with min_idx or max_idx'); end

  % alter branchData
  branchData(branchNr).schnitzNrs = branchData(branchNr).schnitzNrs(min_idx(1):max_idx(end));
  for i = 1:length(p.dataFields)
    field = char(p.dataFields(i));
    branchData(branchNr).(field) = branchData(branchNr).(field)(min_idx(1):max_idx(end));
  end
end

% SHOW SOME OUTPUT
min_branch_length = realmax('double');
max_branch_length = 0;
min_begin_time = realmax('double');
max_begin_time = 0;
min_end_time = realmax('double');
max_end_time = 0;
for branchNr = 1:length(branchData)
  min_branch_length  = min(min_branch_length, length(branchData(branchNr).(timeField)) );
  max_branch_length  = max(max_branch_length, length(branchData(branchNr).(timeField)) );
  min_begin_time  = min(min_begin_time, branchData(branchNr).(timeField)(1) );
  max_begin_time  = max(max_begin_time, branchData(branchNr).(timeField)(1) );
  min_end_time    = min(min_end_time, branchData(branchNr).(timeField)(end) );
  max_end_time    = max(max_end_time, branchData(branchNr).(timeField)(end) );
end
disp(['sameLength = 0']);
disp(['Number of Branches : ' num2str(length(branchData))]);
disp(['Branch time points : from ' num2str(round(min_branch_length)) ' till ' num2str(round(max_branch_length))]);
disp(['Begin time         : from ' num2str(round(min_begin_time)) ' mins till ' num2str(round(max_begin_time)) ' mins']);
disp(['End time           : from ' num2str(round(min_end_time)) ' mins till ' num2str(round(max_end_time)) ' mins']);
disp(['---------------------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% IN CASE SAMELENGTH, REMOVE BRANCHES WITH SHORTER LENGTH
%--------------------------------------------------------------------------
% Check data in case we want sameLength data (remove branches that are
% shorter than timeFrame)
min_begin_time = min ( [branchData.(timeField)] );
max_end_time = max ( [branchData.(timeField)] );

branchDataSameLength = branchData(1);
count = 0;
% loop over branches
for branchNr = 1:length(branchData)
  if branchData(branchNr).(timeField)(1) == min_begin_time & branchData(branchNr).(timeField)(end) == max_end_time
    count = count + 1;
    branchDataSameLength(count) = branchData(branchNr);
  end
end

% SHOW SOME OUTPUT
min_branch_length = realmax('double');
max_branch_length = 0;
min_begin_time = realmax('double');
max_begin_time = 0;
min_end_time = realmax('double');
max_end_time = 0;
for branchNr = 1:length(branchDataSameLength)
  min_branch_length  = min(min_branch_length, length(branchDataSameLength(branchNr).(timeField)) );
  max_branch_length  = max(max_branch_length, length(branchDataSameLength(branchNr).(timeField)) );
  min_begin_time  = min(min_begin_time, branchDataSameLength(branchNr).(timeField)(1) );
  max_begin_time  = max(max_begin_time, branchDataSameLength(branchNr).(timeField)(1) );
  min_end_time    = min(min_end_time, branchDataSameLength(branchNr).(timeField)(end) );
  max_end_time    = max(max_end_time, branchDataSameLength(branchNr).(timeField)(end) );
end
disp(['sameLength = 1']);
disp(['Number of Branches : ' num2str(length(branchDataSameLength))]);
disp(['Branch time points : from ' num2str(round(min_branch_length)) ' till ' num2str(round(max_branch_length))]);
disp(['Begin time         : from ' num2str(round(min_begin_time)) ' mins till ' num2str(round(max_begin_time)) ' mins']);
disp(['End time           : from ' num2str(round(min_end_time)) ' mins till ' num2str(round(max_end_time)) ' mins']);
disp(['---------------------------------------------------------------']);

if p.sameLength
  branchData = branchDataSameLength;
  disp(['Returned sameLength = 1']);
else
  disp(['Returned sameLength = 0']);
end
disp(['---------------------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER BRANCHES AND CALC & SAVE WEIGHT
%--------------------------------------------------------------------------
unique_schnitzNrs = unique([branchData.schnitzNrs]);
unique_timeField  = unique([branchData.(timeField)]);
schnitz_count = zeros(length(unique_schnitzNrs), length(unique_timeField));

% loop over branches
for branchNr = 1:length(branchData)
  for i = 1:length(branchData(branchNr).schnitzNrs)
    idx_schnitzNrs = find(unique_schnitzNrs == branchData(branchNr).schnitzNrs(i));
    idx_timeField  = find(unique_timeField  == branchData(branchNr).(timeField)(i));
    schnitz_count(idx_schnitzNrs, idx_timeField) = schnitz_count(idx_schnitzNrs, idx_timeField) + 1;
  end
end
for branchNr = 1:length(branchData)
  for i = 1:length(branchData(branchNr).schnitzNrs)
    idx_schnitzNrs = find(unique_schnitzNrs == branchData(branchNr).schnitzNrs(i));
    idx_timeField  = find(unique_timeField  == branchData(branchNr).(timeField)(i));
    branchData(branchNr).weight(i) = 1 / schnitz_count(idx_schnitzNrs, idx_timeField);
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% GET MEAN FOR NOISE & NORM FIELDS
%--------------------------------------------------------------------------
datafield_sum = zeros(length(p.dataFields), length(unique_timeField));
datafield_count = zeros(length(p.dataFields), length(unique_timeField));

% loop over branches
for branchNr = 1:length(branchData)
  for age = 1:length(branchData(branchNr).schnitzNrs)
    idx  = find(unique_timeField  == branchData(branchNr).(timeField)(age));
    for i = 1:length(p.dataFields)
      datafield_sum(i,idx) = datafield_sum(i,idx) + branchData(branchNr).(char(p.dataFields(i)))(age)/branchData(branchNr).weight(age);
      datafield_count(i,idx) = datafield_count(i,idx) + 1/branchData(branchNr).weight(age);
    end
  end
end

datafield_mean = datafield_sum ./ datafield_count;

% loop over noiseFields in dataFields
for i = 2:length(p.dataFields)
  field = char(p.dataFields(i));
  noisefield = ['noise_' field];
  normfield  = ['norm_' field];

  % loop over branches
  for branchNr = 1:length(branchData)
    % loop over data
    for age = 1:length(branchData(branchNr).(field))
      idx  = find(unique_timeField  == branchData(branchNr).(timeField)(age));
      branchData(branchNr).(noisefield)(age) = branchData(branchNr).(field)(age) - datafield_mean(i,idx);
      branchData(branchNr).(normfield)(age) = branchData(branchNr).(field)(age) - mean(datafield_mean(i,:));
    end
  end
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% RANKING, for now only works well with sameLength
%--------------------------------------------------------------------------
for i = 2:length(p.dataFields)
  field = char(p.dataFields(i));
  rankfield = ['rank_' field];
  timerankfield = ['timerank_' field];

  for t = 1:length(unique_timeField)
    time = unique_timeField(t);
    data = []; count = 0;

    for branchNr = 1:length(branchData)
      idx_timeField  = find(branchData(branchNr).(timeField) == time);
      if idx_timeField
        if ~isnan(branchData(branchNr).(field)(idx_timeField(1)))
          count = count + 1;
          data(count) = branchData(branchNr).(field)(idx_timeField(1));
        else
          disp([ 'branchNr ' num2str(branchNr) ...
                 ' has schnitzcells(' num2str(branchData(branchNr).schnitzNrs(idx_timeField(1))) ...
                 ').' field ...
                 '(' num2str(idx_timeField(1)) ...
                 ') = NaN (t= ' num2str(time) '). Might give ranking problem']);
        end
      end
    end
    
    ranking = FractionalRankings(data);
    ranking = ranking - mean(ranking);

    count = 0;
    for branchNr = 1:length(branchData)
      idx_timeField  = find(branchData(branchNr).(timeField) == time);
      if idx_timeField
        if ~isnan(branchData(branchNr).(field)(idx_timeField(1)))
          count = count + 1;
          branchData(branchNr).(rankfield)(idx_timeField(1)) = ranking(count);
        else
          branchData(branchNr).(rankfield)(idx_timeField(1)) = NaN;
        end
      end
    end
  end
  
  for branchNr = 1:length(branchData)
    branchData(branchNr).(timerankfield) = FractionalRankings( branchData(branchNr).(field) )';
    branchData(branchNr).(timerankfield) = branchData(branchNr).(timerankfield) - mean( branchData(branchNr).(timerankfield) );
  end
  
end

%--------------------------------------------------------------------------
