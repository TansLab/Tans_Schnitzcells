% function branches = MW_getBranches(p, s, varargin)
%
% MW_getBranches is a modified version of DJK_getBranches, the section how
% a timewindow is selected is updated. (I did this because I encountered
% some issues with the previous version when fitTime(1)>0.)
%
% DJK_getBranches takes schnitzcells and returns branch data of
% particular branches & defined data. 
%
% Data fields that are returned can be set using 'dataFields'. schnitzNrs
% datafield will always be added to branches. 
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
% There is an extra option that only returns branches of a single schnitz
% long (this means each branch contains data from a single cell cycle):
% 'singleSchnitz'. This should be used with 'sameLength'=0
%
% Some data might be used more than once (in different branches). The
% number of branches a particular datapoint is used in, is counted in the
% 'count' field. An additional field 'branchpoints', sometimes used for
% weighing is also added. It counts the number of real branchpoints (real
% meaning, that if the data of the other branch is also taken as
% branches) that occur since the start of the branch.
%
% OUTPUT
% 'branches'          cell structure with branches
%
% REQUIRED ARGUMENTS:
% 'p'
% 's'                 schnitzcells with useForPlot
%
% OPTIONAL ARGUMENTS:
% 'dataFields'        fields to be stored in branches
%                     default: {'Y_time' 'Y6_mean' 'mu_fitNew'}
%                     IMPORTANT NOTES: 
%                       - should have same # values!!!! 
%                         (MW 2015/04 added throwing warning msg in this 
%                         case)
%                         [MW TODO Special case to investigate: what
%                         happens if rates are used? is time then correctly
%                         used and how are fields selected?]
%                       - Put time in first one!
%                         (MW 2015/04 added warning when time not recogn.)
%                       - Some default fields, like mu, contain values that
%                         correspond to the timepoints of when the 
%                         fluor1 images are taken (i.e. this might cause
%                         issues when you are using multiple fluors).
%                     
%
% 'sameLength'=0      does not make branches the same length (default:1)
%
% 'singleSchnitz'=1   makes branch length only 1 schnitz long & not same length (default:0)
%
% 'fitTime'=[10 100]  time between datapoints in branches must lie
%                     (default: all timepoints)
%

function branches = MW_getBranches(p, s, varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'DJK_getBranches';

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
%--------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If not provided, use standard dataFields
if ~existfield(p, 'dataFields')
  p.dataFields = {'Y_time' 'Y6_mean' 'mu_fitNew'};
end

for i = 1:length(p.dataFields)
  field = char(p.dataFields(i));
  if ~existfield(s(1),field)
    disp(['Field ' field ' does not exist. Exiting...!']);
    return;
  end
end

% first dataField should contain time (timeField) (MW 2015/04)
timeField = char(p.dataFields(1));
if isempty(strfind(timeField,'time'))
    warning('The term ''time'' doesn''t occur in the first field''s name; note that this should be a field with time values.');
end

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

%% -------------------------------------------------------------------------
% More input error checking
%--------------------------------------------------------------------------

% Throw error if fields are not same size, since function misbehaves in 
% that case (e.g. schnitzNrs are not correct any more).
% MW 2015/04
numelsDatafields=[];
for i = 1:length(p.dataFields)
  numelsDatafields(end+1) = numel([s.(p.dataFields{i})]); % determine lengths  
end
deltaNumelsDatafields = numelsDatafields-numelsDatafields(1); % this should be zero if they're all equal lengths
% Throw warning if not same size
if any(deltaNumelsDatafields)
    warning(['Given datafields are not same size. This can result in undesired/incorrect behavior! ' ... 
        'Only ignore this warning when you know what you''re doing! (Resuming in 10 seconds.)']); % MW TODO warning not always given??!
    disp('Warning! Given datafields are not same size!'); 
    deltaNumelsDatafields
    pause(10);
end


%% ------------------------------------------------------------------------
% SOME PREPARATION FOR OUTPUT
%--------------------------------------------------------------------------
% remembers whether data for this schnitz has been used already
for i = 1:length(s)
  s(i).schnitzUsed = 0;
end

% branches will contain data
branches = struct;

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
      %empty cell:
      data_temp = cell(1+length(p.dataFields),1);
     
      % now get data for this (starts from schnitzNr) branch from end to beginning
               
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
        % I.e. add N times the current schnitzNr, where N is the length of
        % the timespan of that schnitz's life.
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
               
        %end while:
      end
      
      % invert data_temp 
      for i = 1:length(data_temp)
        data_invers = data_temp{i,:};
        data_temp{i,:} = data_invers(end:-1:1);
      end

      % add to branches
      branches(branchNr).schnitzNrs = data_temp{1,:};
      for i = 1:length(p.dataFields)
        field = char(p.dataFields(i));
        branches(branchNr).(field) = data_temp{i+1,:};
      end

    end
  end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER BRANCHES AND REMOVE DATA OUTSIDE TIME FRAME
%--------------------------------------------------------------------------
% Which fields to process, see below
uniqueFieldsToProcess = unique(p.dataFields); % see below
uniqueFieldsToProcess{end+1} = 'schnitzNrs';
% Loop over branches
for branchNr = 1:length(branches)
    
    % Find appropriate time indices
    timeIndices = (branches(branchNr).(p.dataFields{1})>=p.fitTime(1)) & (branches(branchNr).(p.dataFields{1})<=p.fitTime(2));
    
    % Now apply them to all fields given (but not twice, hence "unique")
    % Apply
    for currentFieldIdx = 1:numel(uniqueFieldsToProcess)
        branches(branchNr).(uniqueFieldsToProcess{currentFieldIdx}) = ...
            branches(branchNr).(uniqueFieldsToProcess{currentFieldIdx})(timeIndices);
    end
        
end


%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% IN CASE SAMELENGTH, REMOVE BRANCHES WITH SHORTER LENGTH
%--------------------------------------------------------------------------

if p.sameLength   
    
    originalSize = numel(branches);
    
    % Throw out empty branches first
    emptyBranchIndices = arrayfun(@(x) isempty(x.(timeField)),branches)
    branches = branches(~emptyBranchIndices);
    
    % get all lengths, begin times, and end times
    branchlengths = arrayfun(@(x) length(x.(timeField)),branches);
    startTimes    = arrayfun(@(x) x.(timeField)(1),branches);
    endTimes      = arrayfun(@(x) x.(timeField)(end),branches);
    % determine min/max values
    longestLength        = max(branchlengths);
    desiredStartTime     = min(startTimes);
    desiredEndTime       = max(endTimes);
        
    % determine which ones are OK
    % length of branch should equal longest 
    condition1 = [branchlengths == longestLength];
    % startime should be desired one (lowest of all)
    condition2 = [startTimes == desiredStartTime];
    % endtime should be desired one (highest of all)
    condition3 = [endTimes == desiredEndTime];
    % and branches to be selected are the ones that satisfy all conditions.
    branchIdxWeWant = condition1 & condition2 & condition3;
    
    % perform the selection
    branches = branches(branchIdxWeWant);
    
    % give user some info
    disp([num2str(numel(branches)) ' branches were selected of total ' num2str(originalSize)]);
    disp(['Selected timewindow is ' num2str(desiredStartTime) ' to ' ...
          num2str(desiredEndTime) ', all with lengths of ' num2str(longestLength)]);
    
  disp(['Returned sameLength = 1']);
else
  disp(['Returned sameLength = 0']);
end
disp(['---------------------------------------------------------------']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER BRANCHES AND CALC & SAVE TIMES DATAPOINT USED (count)
%--------------------------------------------------------------------------
unique_schnitzNrs = unique([branches.schnitzNrs]);
unique_timeField  = unique([branches.(timeField)]);
schnitz_count = zeros(length(unique_schnitzNrs), length(unique_timeField));

% loop over branches
for branchNr = 1:length(branches)
  for i = 1:length(branches(branchNr).schnitzNrs)
    idx_schnitzNrs = find(unique_schnitzNrs == branches(branchNr).schnitzNrs(i));
    idx_timeField  = find(unique_timeField  == branches(branchNr).(timeField)(i));
    schnitz_count(idx_schnitzNrs, idx_timeField) = schnitz_count(idx_schnitzNrs, idx_timeField) + 1;
  end
end
for branchNr = 1:length(branches)
  for i = 1:length(branches(branchNr).schnitzNrs)
    idx_schnitzNrs = find(unique_schnitzNrs == branches(branchNr).schnitzNrs(i));
    idx_timeField  = find(unique_timeField  == branches(branchNr).(timeField)(i));
    branches(branchNr).count(i) = schnitz_count(idx_schnitzNrs, idx_timeField);
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER BRANCHES AND CALC & SAVE BRANCHPOINTS (branchpoints)
%--------------------------------------------------------------------------
% loop over branches
for branchNr = 1:length(branches)
  branches(branchNr).branchpoints(1) = 0;
  for i = 2:length(branches(branchNr).schnitzNrs)
    % in principal same number of branchpoints as timepoint before
    branches(branchNr).branchpoints(i) = branches(branchNr).branchpoints(i-1);

    % check whether there has been a branch, then add 1
    % first check whether sister schnitz exists in any of the branches (to
    % check whether it is a real branch point)
    idx = find(s(branches(branchNr).schnitzNrs(i)).S == unique_schnitzNrs);
    if branches(branchNr).schnitzNrs(i) ~= branches(branchNr).schnitzNrs(i-1) & ...
       length(idx)>0
      branches(branchNr).branchpoints(i) = branches(branchNr).branchpoints(i) + 1;
    end
  end
end
%--------------------------------------------------------------------------

