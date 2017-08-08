function [dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = MW_getdelayedscatter(p, branches, fieldX, fieldY, allowedRedundancy, varargin);
% function [dataPairsPerTau, iTausCalculated, originColorPerTau, correlationsPerTau] = MW_getdelayedscatter(p, branches, fieldX, fieldY, allowedRedundancy, varargin);
%
% [MW todo: update function declaration]
% This function is based on DJK_getcrossCor.m.
% MW 2015/04
%
% [MW todo: describe input/output of this function!]
% INPUT
% This function requires non-normalized data, i.e. "raw" branches, as
% opposed to the (cross-)corr functions. <-> Does it? Maybe not??!
%
% - allowedRedundancy   This parameter sets how many times datapoints are
%                       allowed to be used again. (Since a time delay is
%                       involved, and branches are redundant, redundancy
%                       cannot be avoided.)
%                       Not sure how necessary this parameter is, since we
%                       also require that one of the points in the
%                       datapair is unique.
%
% OUTPUT
% - dataPairsPerTau     data pairs per tau value, tau values are listed in
%                       iTausCalculated. 
% - originColorPerTau   originColorPerTau gives a color which relates to 
%                       which branch (=lineage) was used for this datapair 
%                       (indexing same as dataPairsPerTau).
%
% Optional arguments:
% tauIndices    Custom range of taus to consider.


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 5; functionName = 'MW_getdelayedscatter';

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

% =========================================================================

% Get total # branches
numelBranches = numel(branches);

% Create color map to identify branches
lineageColormap = colormap(lines(numelBranches));

% loop over different time-lags tau
% MW TODO: Note that this will in some cases produce a rather massive
% dataset, so maybe it is more convenient to give a restraint to t_max
% manually, or better a manual range for tauIndices by the user.
dataPairsPerTau = {}; originColorPerTau = {};
iTausCalculated = []; 
correlationsPerTau = [];
% Loop over iTau, i.e. the delay expressed in terms of the index.
for iTau = p.tauIndices
    
  % Keep track of which schnitzes have been used, reset this per tau value
  schnitzUseRegister = zeros(1,max([branches.schnitzNrs]));  
    
  % Make list of iTau values treated (in principle this should be known by
  % the user already).
  iTausCalculated(end+1) = iTau;
  % Clear list of pairs from previous iTau value
  pairCollectionForTau = []; originCollectionForTau = [];
  
  % Loop over branches, adding data pairs X(t),Y(t+tau) only when they are
  % not redundant (or below a redundancy treshold).
  for br = 1:numelBranches
      
    % Load data in more conveniently named parameters
    X = branches(br).(fieldX);
    Y = branches(br).(fieldY);
    schnitzesInBranch = branches(br).schnitzNrs;    

    %length of this particular branch
    N = length(X);

    % Collect the datapairs
    % Loop over entries in branche in reverse such that least redundant 
    % part of the branche is analyzed first, take lag window in account
    endOfDomain = N; 
    if iTau > 0, endOfDomain = endOfDomain-iTau; end % handle positive tau
    startOfDomain = 1;
    if iTau < 0, startOfDomain = startOfDomain-iTau; end % handle negative tau
    schnitzUsedPastLoop1 = 0; schnitzUsedPastLoop2 = 0;
    for iTime = endOfDomain:-1:startOfDomain % reverse loop
        % note the loop is reversed since we want to back-track the
        % lineage; This changes the chronology, but I don't think that
        % matters anywhere.

        % which schnitzes are we going to use?
        schnitzUsed1 = schnitzesInBranch(iTime);
        schnitzUsed2 = schnitzesInBranch(iTime+iTau);

        % Is this allowed?
        % Check only if we are entering a domain belonging to another
        % schnitz. Since we do want multiple contributions per schnitz.
        if (schnitzUsedPastLoop1~=schnitzUsed1) | ...
           (schnitzUsedPastLoop2~=schnitzUsed2)
       
            % Check 1: one of the contributions should be unique
            if schnitzUseRegister(schnitzUsed1)>1 & ...
               schnitzUseRegister(schnitzUsed2)>1
                % if both are already used, this would be a duplicate
                % datapoint, we don't want that.
                % Note that we use ">1" since a calculation over one
                % branch will hit a schnitz twice.
                break;
            end
            
            % Check 2: neither may be used too much
            if schnitzUseRegister(schnitzUsed1)>allowedRedundancy | ...
               schnitzUseRegister(schnitzUsed2)>allowedRedundancy  

                % if (since we're going backwards in the branches) we have
                % reached a point we're the redundancy is already too high,
                % we should stop running over this branch.
                break;
            end

            % Increase in a register how many times both of these schnitzes 
            % have already been used.
            schnitzUseRegister(schnitzUsed1) = schnitzUseRegister(schnitzUsed1)+1;
            schnitzUseRegister(schnitzUsed2) = schnitzUseRegister(schnitzUsed2)+1;
            
            
        end

        % How redundant is this datapoint?
        % MW TODO!
        % redundancyDatapoint = branches(br).count(XX) + branches(br).count(XX);
        
        % Color code for datapoint to know from which lineage it came
        currentOriginColor = lineageColormap(br,:);
        
        % collect a pair
        currentPair = [X(iTime)  Y(iTime+iTau)];

        % Save this pair
        pairCollectionForTau = [pairCollectionForTau; currentPair];
        % And save its color code for the origin
        originCollectionForTau = [originCollectionForTau; currentOriginColor];
        
        % Remember we used these schnitzes
        schnitzUsedPastLoop1 = schnitzUsed1;
        schnitzUsedPastLoop2 = schnitzUsed2;
        
    end % end time loop

  end % end branch loop

  % in case chronology does matter, one could fix this by reversing the
  % order (untested code)
  %pairCollectionForTau=fliplr(pairCollectionForTau);
  %originCollectionForTau=fliplr(originCollectionForTau);
  
  % Calculate covariance for this tau value
  correlationThisTau = corr(pairCollectionForTau(:,1),pairCollectionForTau(:,2));
  
  % Save data for this tau value
  dataPairsPerTau{end+1} = pairCollectionForTau;
  originColorPerTau{end+1} = originCollectionForTau;
  correlationsPerTau(end+1) = correlationThisTau; 
  
  disp(['Now at ' num2str(iTau) ', going to ' max(num2str(p.tauIndices)) '.']);
  
end % end delay loop (tau)
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
