function [dataPairsPerTau, iTausCalculated] = MW_getdelayedscatter(p, branches, fieldX, fieldY, varargin);
% [MW todo: add function declaration]
% This function is based on DJK_getcrossCor.m.
% MW 2015/04
%
% [MW todo: describe input/output of this function!]
% INPUT
% This function requires non-normalized data, i.e. "raw" branches, as
% opposed to the (cross-)corr functions.
% 
% Optional arguments:
% tauIndices    Custom range of taus to consider.


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_getCrossCov';

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
% Set p struct fields that are not set to default values.
% --------------------------------------------------------------------------
if ~existfield(p,'timeField')
  disp('p.timeField not set, reverting to default ''Y_time''.');
  p.timeField = 'Y_time';
end

% check whether fields exist
if ~existfield(branches(1),fieldX)
  disp(['Field ' fieldX ' does not exist. Exiting...!']); return;
end
if ~existfield(branches(1),fieldY)
  disp(['Field ' fieldY ' does not exist. Exiting...!']); return;
end
if ~existfield(branches(1), p.timeField)
  disp(['Field ' p.timeField ' does not exist (defined in p.timeField). Exiting...!']); return;
end

% if too little data, exit
if (length(branches(1).(fieldX)) < 2), error('Not enough data'); end

% Some important settings
if ~existfield(p,'spacingError')
  p.spacingError = 0.05;
end

% Defuaulf value for field tauIndices is set below.

% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Check whether branches have same length
% --------------------------------------------------------------------------
maxBranchLength = func_check_branchLength(branches);
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Check sampling interval for even spacing
% --------------------------------------------------------------------------
interval = func_check_spacingError(branches(1).(p.timeField), p.spacingError);
% time for crosscorrelation
for br = 1:length(branches)
    % Calculate a time vector from [-t_max:dt:t_max], note that this vector
    % only corresponds to the real vector within the margin of the spacing 
    % error if there were no warnings.
    branches(br).corr_time = interval * [-maxBranchLength+1:maxBranchLength-1];    
end
% --------------------------------------------------------------------------

% Delays, in terms of index, to be taken into account
if ~isfield(p,'tauIndices')
    p.tauIndices = [-maxBranchLength+1:maxBranchLength-1];
end

% loop over different time-lags tau
% MW TODO: Note that this will in some cases produce a rather massive
% dataset, so maybe it is more convenient to give a restraint to t_max
% manually, or better a manual range for tauIndices by the user.
dataPairsPerTau = struct;
iTausCalculated = [];
schnitzUseRegister = []
% Loop over iTau, i.e. the delay expressed in terms of the index.
for iTau = p.tauIndices
  % Make list of iTau values treated (in principle this should be known by
  % the user already).
  iTausCalculated(end+1) = iTau;
  % Clear list of pairs from previous iTau value
  pairCollectionForTau = [];
  
  % Loop over branches, adding data pairs X(t),Y(t+tau) only when they are
  % not redundant (or below a redundancy treshold).
  for br = 1:length(branches)      
      
    % Load data in more conveniently named parameters
    X = branches(br).(fieldX);
    Y = branches(br).(fieldY);
    schnitzesInBranch = branches(br).schnitzNrs;

    %length of this particular branch
    N = length(X);

    % Collect the datapairs
    % Loop over entries in branche in reverse such that least redundant 
    % part of the branche is analyzed first
    for iTime = N:-1:1
        
        % MW TODO: only when schnitzes not used too many times before
        if 1 % MW TODO MAKE THIS IF STATEMENT! 
            % collect a pair
            currentPair = [X(iTime)  Y(iTime+iTau)]

            % Save this pair         
            pairCollectionForTau = [pairCollectionForTau; currentPair];
                
            % Increase in a register how many times both of these schnitzes 
            % have already been used.
            schnitzUsed1 = schnitzesInBranch(iTime)
            schnitzUsed2 = schnitzesInBranch(iTime+iTau)
        else
            % break the loop, since we're entering redundant territory
            break;
        end
        
    end
    
    %{
    % in case lag is positive
    if iTau >= 0
        
        % WHERE MAGIC SHOULD HAPPEN (1)
        X(iTau) .* Y(1+iTau:N)
        
      %branches(br).(targetField)(N+iTau) = sum( X(1:N-iTau) .* Y(1+iTau:N) );
      %crossCov_composite.Y(N+iTau) = crossCov_composite.Y(N+iTau) + sum( X(1:N-iTau) .* Y(1+iTau:N) .* W );
    
    % in case lag is negative
    else
        
        % WHERE MAGIC SHOULD HAPPEN (2)        
        
        
      %branches(br).(targetField)(N+iTau) = sum( Y(1:N+iTau) .* X(1-iTau:N) );
      %crossCov_composite.Y(N+iTau) = crossCov_composite.Y(N+iTau) + sum( Y(1:N+iTau) .* X(1-iTau:N) .* W );
    end
    %}
    
  end % end branch loop

end % end delay loop
% -------------------------------------------------------------------------

% Now final data should be produced, output struct is three dimensional:
% - output(iTau) gives array with X(t),Y(t+tau) datapairs.
% - output(iTau,j) gives pair j
% - output(iTau,j,k) gives pair j, observable k (k=1 -> X; k=2 -> Y).

% --------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = func_check_spacingError(timeData, spacingError), ...
% Check whether timefield values are evenly spaced. 
% (Since tau is calculated as dt*i, where i is the index of the field of
% interest. Hence, dt, the spacing between the entries, should be 
% constant.)

timeData = timeData-timeData(1);
fit = polyfit([1:length(timeData)],timeData,1);
% disp([' * Time interval is ' num2str(fit(1)) ' (mins)']); %blubb
% NW2012-10

fittedTimeData = fit(2) + fit(1)*[1:length(timeData)];
interval_error = 0;
for i = 1:length(timeData)
  if abs(timeData(i) - fittedTimeData(i)) > spacingError*fit(1)
    if ~interval_error, warning([' * Watch out! timepoints might not be evenly spaced! See graph.']); end
    interval_error = 1;
    disp(['   -> at i=' num2str(i) ' : measured time = ' num2str(timeData(i)) ' & fitted time = ' num2str(fittedTimeData(i))]);
  end
end
if interval_error
  figure; plot(timeData,'bo'); hold on; plot(fittedTimeData,'k-'); 
end
out = fit(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether branches have same length
function max_br_length = func_check_branchLength(branches), ...
max_br_length = 0;
min_br_length  = realmax('double');
for br = 1:length(branches)
  max_br_length = max(max_br_length, length(branches(br).schnitzNrs));
  min_br_length = min(min_br_length, length(branches(br).schnitzNrs));  
end

if (max_br_length == min_br_length)
  %disp([' * All branches have same length: ' num2str(min_br_length) '
  %timepoints.']); blubb NW2012-010 just commented out
else
  disp([' * Branches have different lengths: from ' num2str(min_br_length) ' to ' num2str(max_br_length) ' timepoints. Watch out! Do not know whether composite is correct!']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
