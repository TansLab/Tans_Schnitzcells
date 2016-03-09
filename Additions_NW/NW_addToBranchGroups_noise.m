function branchGroupsNoise = NW_addToBranchGroups_noise(p, branchGroups, varargin)
%
% *********************************************************************
% Follows strongly Daan's version DJK_addToBranches_noise.
% Difference: Here, noise is added after the Branches are divided into
% (typically 4) BranchGroups. The mean at each time point, which is
% subtracted to get the 'noise' term, is taken seperately for each
% individual branch (total independency of branches).
% Daan version: mean is taken from all 4 branches together
% On long time scales this should in theory end up to the same ... we'll
% see...
%
% *********************************************************************
%
% Of each field except the timefield (dataFields(1)) a noise and norm
% version are added.
%
% noise: normalized by subtracting mean for this timepoint off all branches. 
%       New: The mean is taken seperately for each individual branch
%
% norm: normalized by subtracting the mean of the branch data for all timepoints.  
%       New: temporal mean of each branchGroup seperately
%
% OUTPUT
% 'branchGroupsNoise'    cell structure with brancheGroups, with noise/norm added
%
% REQUIRED ARGUMENTS:
% 'branchGroups'         cell structure with branchGroups
%
% OPTIONAL ARGUMENTS:
% 'dataFields'           fields to be stored in branches
%                        default: {'Y_time' 'Y6_mean' 'mu_fitNew'}
%


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'NW_addToBranchGroups_noise';

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
% If not provided, use standard dataFields
if ~existfield(p, 'dataFields')
  p.dataFields = {'Y_time' 'Y6_mean' 'mu_fitNew'};
end

for i = 1:length(p.dataFields)
  field = char(p.dataFields(i));
  if ~existfield(branchGroups(1).branches(1),field) % new check when branchGroup instead of branches are input (NW 2012/03)
    disp(['Field ' field ' does not exist. Exiting...!']);
    return;
  end
end

% first dataField should contain time (timeField)
timeField = char(p.dataFields(1));

% loop over all branchGroups and get unique time fields (NW 2012/03)
clear unique_timeField
for i=1:length(branchGroups)
    unique_timeField(i,:)  = unique([branchGroups(i).branches.(timeField)]);
end
unique_timeField=unique(unique_timeField);
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% LOOP OVER ALL BRANCH GROUPS AND ADD NORM AND NOISE SEPERATELY FOR EACH
% BRANCHGROUP
%--------------------------------------------------------------------------
for run=1:length(branchGroups)
    clear branches datafield_sum datafield_count datafield_mean mybranches;
    mybranches=branchGroups(run).branches;

    %--------------------------------------------------------------------------
    % GET MEAN FOR NOISE & NORM FIELDS
    %--------------------------------------------------------------------------
    datafield_sum = zeros(length(p.dataFields), length(unique_timeField));
    datafield_count = zeros(length(p.dataFields), length(unique_timeField));

    % loop over branches
    for branchNr = 1:length(mybranches)
      for age = 1:length(mybranches(branchNr).schnitzNrs)

        idx  = find(unique_timeField  == mybranches(branchNr).(timeField)(age)); %single index (typically runs from 1 to length of branch= # time points with increment=1, that is age=idx) (comment NW)
        for i = 1:length(p.dataFields)
          datafield_sum(i,idx) = datafield_sum(i,idx) + mybranches(branchNr).(char(p.dataFields(i)))(age) / mybranches(branchNr).count(age); % DJK 091125 last / used to be *, but I now think this is wrong [Yes (NW)]
          datafield_count(i,idx) = datafield_count(i,idx) + 1/mybranches(branchNr).count(age); % DJK 091125 last + 1/ used to be + , but I now think this is wrong [Yes (NW).  I think datafield_count(idx) states total number of schnitzes/cells at time idx)
        end
      end
    end

    datafield_mean = datafield_sum ./ datafield_count;%[comment NW: I think this is the 'normal' mean: Sum(YFP-intensities of all cells)/(# cells in this frames). I'm pretty sure that the weighing with 'count' is cancelled out]

    % loop over noiseFields in dataFields
    for i = 2:length(p.dataFields)
      field = char(p.dataFields(i));
      noisefield = ['noise_' field];
      normfield  = ['norm_' field];

      % loop over branches
      for branchNr = 1:length(mybranches)
        % loop over data
        for age = 1:length(mybranches(branchNr).(field))
          idx  = find(unique_timeField  == mybranches(branchNr).(timeField)(age));
          mybranches(branchNr).(noisefield)(age) = mybranches(branchNr).(field)(age) - datafield_mean(i,idx);
          mybranches(branchNr).(normfield)(age) = mybranches(branchNr).(field)(age) - mean(datafield_mean(i,:));
        end
      end
    end
    %--------------------------------------------------------------------------
    
    % write branches back into branchGroup format
    branchGroups(run).branches=mybranches;
    branchGroupsNoise=branchGroups; % not really necessary. change output variable (NW)
    %--------------------------------------------------------------------------
end
