% function branches = DJK_getBranches(p, s, varargin)
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

function branches = DJK_getBranches(p, s, varargin)

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


%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
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
% remove branches which don't overlap at all with fitTime. Can affect
% weighing! NW2012-09-12
removeBranchNrs=[];
for br=1:length(branches)
    if max(branches(br).(timeField))<p.fitTime(1);
        removeBranchNrs=[removeBranchNrs, br];
        disp('')
        disp(['Will delete short branch ' num2str(br) '. Can affect weighing.']);
    end
end
branches(removeBranchNrs)=[];
% end NW


% DETERMINE TIME FRAME FOR ALL BRANCHES
% Throw all recored times on a pile
all_time_data = sort([branches.(timeField)]);
% Find which is closest to the given boundary
idx = find(all_time_data >= p.fitTime(1));
% Store that one
timeFrame(1) = all_time_data(idx(1));
% Idem for other boundary
idx = find(all_time_data <= p.fitTime(2));
timeFrame(2) = all_time_data(idx(end));

disp(['---------------------------------------------------------------']);
disp(['fitTime            : ' num2str(round(p.fitTime(1))) ' mins till ' num2str(round(p.fitTime(2))) ' mins']);
disp(['useForPlot time    : ' num2str(round(min(all_time_data))) ' mins till ' num2str(round(max(all_time_data))) ' mins']);
disp(['timeFrame          : ' num2str(round(timeFrame(1))) ' mins till ' num2str(round(timeFrame(2))) ' mins']);
disp(['---------------------------------------------------------------']);

%branches=branches(2:end); % blubb NW2012-05 (manually delete branches that collide with fitTime (errormessage in loop below)
% loop over branches
for branchNr = 1:length(branches)
  
  min_idx = find(branches(branchNr).(timeField)>=timeFrame(1));
  max_idx = find(branches(branchNr).(timeField)<=timeFrame(2));
  
  if length(min_idx)==0 | length(max_idx)==0 
      disp('error with min_idx or max_idx'); 
  end
  %branchNr %debug
  
  startTimeWindowIdx=min_idx(1);
  endTimeWindowIdx=max_idx(end);  
  
  % alter branches
  branches(branchNr).schnitzNrs = branches(branchNr).schnitzNrs(startTimeWindowIdx:endTimeWindowIdx);
  % This can lead to e.g. wrong time fields if a timefield with more
  % entries (e.g. R_time) is cut to a time field with less entries (e.g.
  % dR5_time)
  % ***
  % and: if errormessage: 2 options (dependeing on error):
  % - if one branch is too short, delete it out of the "branches"right
  % before this loop
  % - if fitTime(2) is too high (larger than framerange), rates will be one
  % entry short compared to concentrations. lower fitTime(2) a bit. (NW
  % 2012-06-15)
  for i = 1:length(p.dataFields)
    field = p.dataFields{i};
    
      % debug
      if endTimeWindowIdx>numel(branches(branchNr).(field))
          warning('Issue detected! Only proceed with caution! Check whether this is correct!');
          
         % Debugging:     
         disp(['branchNr = ' num2str(branchNr)]); 
         disp(['field = ' field]); 
         disp(['min_idx(1) = ' num2str(min_idx(1))]); 
         disp(['max_idx(end) = ' num2str(max_idx(end))]); 
         disp(['length(branches(branchNr).(field)) = ' num2str(length(branches(branchNr).(field)))]); 
         disp(['----------------------------------------------------------']);      
         
         disp('Setting endTimeWindowIdx to numel(branches(branchNr).(field))..');
         %pause(4);
         
         % This type of behavior (see also comment NW above) happens when
         % fields have incosistent sizes).. 
         % If this is the case for a load of branches, your data will get
         % screwed..
         endTimeWindowIdx = numel(branches(branchNr).(field));
      end
          
    branches(branchNr).(field) = branches(branchNr).(field)(startTimeWindowIdx:endTimeWindowIdx);
  end
end

% SHOW SOME OUTPUT
min_branch_length = realmax('double');
max_branch_length = 0;
min_begin_time = realmax('double');
max_begin_time = 0;
min_end_time = realmax('double');
max_end_time = 0;
for branchNr = 1:length(branches)
  min_branch_length  = min(min_branch_length, length(branches(branchNr).(timeField)) );
  max_branch_length  = max(max_branch_length, length(branches(branchNr).(timeField)) );
  min_begin_time  = min(min_begin_time, branches(branchNr).(timeField)(1) );
  max_begin_time  = max(max_begin_time, branches(branchNr).(timeField)(1) );
  min_end_time    = min(min_end_time, branches(branchNr).(timeField)(end) );
  max_end_time    = max(max_end_time, branches(branchNr).(timeField)(end) );
end
disp(['sameLength = 0']);
disp(['Number of Branches : ' num2str(length(branches))]);
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
min_begin_time = min ( [branches.(timeField)] );
max_end_time = max ( [branches.(timeField)] );

branchesSameLength = branches(1);
count = 0;
% loop over branches
for branchNr = 1:length(branches)
  if branches(branchNr).(timeField)(1) == min_begin_time & branches(branchNr).(timeField)(end) == max_end_time
    count = count + 1;
    branchesSameLength(count) = branches(branchNr);
  end
end

% SHOW SOME OUTPUT
min_branch_length = realmax('double');
max_branch_length = 0;
min_begin_time = realmax('double');
max_begin_time = 0;
min_end_time = realmax('double');
max_end_time = 0;
for branchNr = 1:length(branchesSameLength)
  min_branch_length  = min(min_branch_length, length(branchesSameLength(branchNr).(timeField)) );
  max_branch_length  = max(max_branch_length, length(branchesSameLength(branchNr).(timeField)) );
  min_begin_time  = min(min_begin_time, branchesSameLength(branchNr).(timeField)(1) );
  max_begin_time  = max(max_begin_time, branchesSameLength(branchNr).(timeField)(1) );
  min_end_time    = min(min_end_time, branchesSameLength(branchNr).(timeField)(end) );
  max_end_time    = max(max_end_time, branchesSameLength(branchNr).(timeField)(end) );
end
disp(['sameLength = 1']);
disp(['Number of Branches : ' num2str(length(branchesSameLength))]);
disp(['Branch time points : from ' num2str(round(min_branch_length)) ' till ' num2str(round(max_branch_length))]);
disp(['Begin time         : from ' num2str(round(min_begin_time)) ' mins till ' num2str(round(max_begin_time)) ' mins']);
disp(['End time           : from ' num2str(round(min_end_time)) ' mins till ' num2str(round(max_end_time)) ' mins']);
disp(['---------------------------------------------------------------']);

if p.sameLength
  branches = branchesSameLength;
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

