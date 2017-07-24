function [branches, crossCov_composite] = DJK_getCrossCov(p, branches, fieldX, fieldY, varargin);
% [branches, crossCov_composite] = DJK_getCrossCov(p, branches, fieldX, fieldY, varargin)
%
% DJK_getCrossCov returns the cross covariance of a set of branches,
% together with the composite cross covariance. Different settings can be
% set, such as: 
%  * not biased or biased (Elowitz vs Simpson)
%  * additional scaling of branch (subtracting mean of each branch)
%  * for the composite different weighing methods
%
% The auto cross covariance can be obtained, by providing the dataField
% twice (both as fieldX and fieldY). Then only [0 ->] is returned
%
% A check is performed whether the time points are evenly spaced. using 
% 'spacingError' you can set the maximum allowed offset from even spacing,
% before a warning is shown.
%
% % uses Formula:   r * interval = tau
% %                 Sxy(r) = n=1_SUM_N-r ( X(t) * Y(t+r) ) / (N-|r|)
%
% biased:         in Sxy(r) don't divide by (N-|r|) -> Nature439_608 way
%                 for larger delay times the autocorrelation will always decrease
%
% OUTPUT
% 'branches'            branches structure with calculated cross covariance for each branch added
% 'crossCov_composite'  struct with Y = calculated composite cross covariance for branches
%                                   X = corr time
%
% REQUIRED ARGUMENTS:
% 'p'
% 'branches'
% 'fieldX'        field X of which crosscorrelation should be calculated
% 'fieldY'        field Y of which crosscorrelation should be calculated
%
% OPTIONAL ARGUMENTS:
% 'timeField'     field used as time and stored in branches
%                 default: 'Y_time'
% 'bias'=1        do not adjust for less data at larger delay times (default: 0)
% 'weighing'  =0  performs no weighing for composite 
%             =1  performs standard weighing (like Elowitz) (default=1)
%             =2  performs 3/4 weighing
%             =3  performs Daan weighing
% 'spacingError'  error allowed in spacing (in fraction of average spacing)
%                 default: 0.05
% 'extraNorm'=1   perform an extra normalization where the mean of each
%                 branch is subtracted (default: 0)
% 'rzeroonly'=1    used for crossCor when only r=0 is needed
% 'sameLength'    to signal that function should handle branches of uneven
%                 length

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
% overwrite any schnitzcells parameters/defaults given optional fields/values
% --------------------------------------------------------------------------
%%
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

% settings
if ~existfield(p,'bias')
  p.bias = 0;
end
if ~existfield(p,'weighing')
  p.weighing = 1;
end
if ~existfield(p,'spacingError')
  p.spacingError = 0.05;
end
if ~existfield(p,'extraNorm')
  p.extraNorm = 0;
end
if ~existfield(p,'rzeroonly')
  p.rzeroonly = 0;
end
% --------------------------------------------------------------------------
%%
if ~p.rzeroonly
  % --------------------------------------------------------------------------
  % Check whether branches have same length
  % --------------------------------------------------------------------------  
  maxBranchLength=max(arrayfun(@(x) numel(x.schnitzNrs),branches));
  % %maxBranchLength = func_check_branchLength(branches); % old way - sub function
  % --------------------------------------------------------------------------

    %%
  % --------------------------------------------------------------------------
  % Check sampling interval for even spacing
  % --------------------------------------------------------------------------
  allTimesUnique = unique([branches(:).(p.timeField)]);
  interval = func_check_spacingError(allTimesUnique, p.spacingError);

  %%
  % time for crosscorrelation
  for br = 1:length(branches)
    branches(br).corr_time = interval * [-maxBranchLength+1:maxBranchLength-1];
  end
  % --------------------------------------------------------------------------

  rs = [-maxBranchLength+1:maxBranchLength-1];

else
  maxBranchLength = length(branches(1).schnitzNrs);
  rs = [0];
  interval = 0;
end
%%
% --------------------------------------------------------------------------
% Calc of correlation for each branch
% --------------------------------------------------------------------------
% crosscorrelation for each independent branch will be saved in branches in this field
targetField = ['crossCov_' fieldX '_' fieldY];
%blubb % force truncation to 63 characters to suppress warning NW
%2012-05-16. AKWARD!!!!
if length(targetField)>63
    length(targetField);
    targetField=targetField(1:63);
    length(targetField);
end
% blubb ende

% prepare composite
crossCov_composite = struct;
crossCov_composite.Y = zeros(1,length([-maxBranchLength+1:maxBranchLength-1]));
crossCov_composite.X = interval * [-maxBranchLength+1:maxBranchLength-1];

%% loop over different time-lags r

%total_weight = nan(1,numel(rs));

for r = rs

  total_weight = 0; % calc total weight for this r

  % loop over branches
  for br = 1:length(branches)           

    X = branches(br).(fieldX);
    Y = branches(br).(fieldY);

    % in case extra normalization, subtract mean of branch
    if p.extraNorm
      X = X - mean(X);
      Y = Y - mean(Y);
    end

    %length of this particular branch
    N = length(X);

    % determine weighing
    clear W;
    switch p.weighing
      case 0 % performs no weighing for composite 
        W = ones(1,N-abs(r));
      case 1 % weighing 1 performs standard weighing (like Elowitz) (default)
        W = 1 ./ branches(br).count(1+abs(r):end); 
            % NOTE MW: Here, only the count of one point contributing to
            % pair is used, because this is equivalent to subtracting
            % non-unique pairs (a pair of points is already redundant
            % if only one point being redundant).
      case 2 % weighing 2 performs 3/4 weighing
        W = 1 ./ branches(br).count(1+abs(r):end); 
        for i = 1:N-abs(r)
          if branches(br).count(i) > 1
            W(i) = W(i)*0.75;
          end
        end
      case 3 % weighing 3 performs Daan weighing
        for i = 1:N-abs(r)
          B = branches(br).branchpoints(i+abs(r)) - branches(br).branchpoints(i);
          W(i) = (1/(2*branches(br).count(i+abs(r)))) * (1 + 2^B);
        end
      otherwise
        disp('Unknown weighing method! Using default');
        W = 1 ./ branches(br).count(1+abs(r):end); 
    end

    total_weight = total_weight + sum( W );

    % in case r is positive
    if r >= 0
      branches(br).(targetField)(N+r) = sum( X(1:N-r) .* Y(1+r:N) );
      crossCov_composite.Y(N+r) = crossCov_composite.Y(N+r) + sum( X(1:N-r) .* Y(1+r:N) .* W );
    % in case r is negative
    else
      branches(br).(targetField)(N+r) = sum( Y(1:N+r) .* X(1-r:N) );
      crossCov_composite.Y(N+r) = crossCov_composite.Y(N+r) + sum( Y(1:N+r) .* X(1-r:N) .* W );
    end

    % if unbiased divide by (N-r) -> Elowitz is unbiased, Simpson is biased
    if ~p.bias
      branches(br).(targetField)(N+r) = branches(br).(targetField)(N+r) / (N-abs(r));
    end

  end

  % divide by total weight
  crossCov_composite.Y(maxBranchLength+r) = crossCov_composite.Y(maxBranchLength+r) / total_weight;

  % if unbiased divide by (N-r) -> Elowitz is unbiased, Simpson is biased
  if p.bias
    crossCov_composite.Y(maxBranchLength+r) = crossCov_composite.Y(maxBranchLength+r) * (N-abs(r));
  end
end
% -------------------------------------------------------------------------

%%
% --------------------------------------------------------------------------
% In case of autocorrelation, only return [0 ->]
% --------------------------------------------------------------------------
% check whether we're dealing with an autocorrelation
% if p.rzeroonly, already doing only 0
if strcmp(fieldX, fieldY) & ~p.rzeroonly
  
  % Find location of zero on x-axis (tau)
  idx = find(crossCov_composite.X==0); %MW 2014/07/30 neater code
  
  % Select only values from X=0 to X(end).
  crossCov_composite.X = crossCov_composite.X(idx:end); % MW 2014/07/30 neater code
  crossCov_composite.Y = crossCov_composite.Y(idx:end); % MW 2014/07/30 neater code

  % now same as above for all branches
  for br = 1:length(branches)
      
    % Find location of zero on x-axis (tau)
    idx = find(branches(br).corr_time==0); % MW 2014/07/30 neater code
    
    % Select only values from X=0 to X(end).
    branches(br).corr_time = branches(br).corr_time(idx:end); % MW 2014/07/30 neater code
    branches(br).(targetField) = branches(br).(targetField)(idx:end); % MW 2014/07/30 neater code
       
  end 

    disp('Autocorr detected.');
  
end
% --------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether constant sampling in timefield
function interval = func_check_spacingError(timeData, spacingError), ...
%timeData=allTimesUnique, spacingError= p.spacingError
%%
timeData = timeData-timeData(1);
fit = polyfit([1:length(timeData)],timeData,1);
% disp([' * Time interval is ' num2str(fit(1)) ' (mins)']); %blubb
% NW2012-10

fittedTimeData = fit(2) + fit(1)*[1:length(timeData)];
interval_error = 0;
for i = 1:length(timeData)
  if abs(timeData(i) - fittedTimeData(i)) > spacingError*fit(1)
    if ~interval_error, disp([' * Watch out! timepoints might not be evenly spaced! See graph [disabled in code]']); end
    interval_error = 1;
    disp(['   -> at i=' num2str(i) ' : measured time = ' num2str(timeData(i)) ' & fitted time = ' num2str(fittedTimeData(i))]);
  end
end
if interval_error
  %figure; plot(timeData,'bo'); hold on; plot(fittedTimeData,'k-'); %enable
  %for investigation. otherwise slows matlab down
end
interval = fit(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% check whether branches have same length
function maxBranchLength = func_check_branchLength(branches), ...
%%
maxBranchLength = 0;
min_br_length  = realmax('double');
for br = 1:length(branches)
  maxBranchLength = max(maxBranchLength, length(branches(br).schnitzNrs));
  min_br_length = min(min_br_length, length(branches(br).schnitzNrs));  
end

if (maxBranchLength == min_br_length)
  %disp([' * All branches have same length: ' num2str(min_br_length) '
  %timepoints.']); blubb NW2012-010 just commented out
else
  disp([' * Branches have different lengths: from ' num2str(min_br_length) ' to ' num2str(maxBranchLength) ' timepoints. Watch out! Do not know whether composite is correct!']);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
