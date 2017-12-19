function [branches, crossCor_composite] = MW_getCrossCor(p, branches, fieldX, fieldY, varargin);
% function [branches, crossCor_composite] = MW_getCrossCor(p, branches, fieldX, fieldY, varargin);
%
% DJK_getCrossCor returns the cross correlation of a set of branches,
% together with the composite cross correlation. Different settings can be
% set, such as: 
%  * not biased or biased (Elowitz vs Simpson)
%  * additional scaling of branch (subtracting mean of each branch)
%  * for the composite different weighing methods
%
% The auto cross correlation can be obtained, by providing the dataField
% twice (both as fieldX and fieldY)
%
% A check is performed whether the time points are evenly spaced. using 
% 'spacingError' you can set the maximum allowed offset from even spacing,
% before a warning is shown.
%
% biased:         in Sxy(r) don't divide by (N-|r|) ->  way
%                 for larger delay times the autocorrelation will always decrease
%
% OUTPUT
% 'branches'            branches structure with calculated cross covariance for each branch added
% 'crossCor_composite'  struct with Y = calculated composite cross correlation for branches
%                                   X = corr time
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
%                 branch is subtracted (default: 0). Note that
%                 normalization of branches should in principle already be
%                 tackled by substracting the whole-colony average at
%                 that timepoint. (One should use only noise_.. fields to
%                 perform crosscorrelations on when using at branches.)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'MW_getCrossCor';

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


% --------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
% --------------------------------------------------------------------------
if ~existfield(p,'timeField')
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
  disp(['Field ' p.timeField ' does not exist. Exiting...!']); return;
end

% if too little data, exit
if (length(branches(1).(fieldX)) < 2), warning('Not enough data'); end

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
% --------------------------------------------------------------------------


% --------------------------------------------------------------------------
% Calc of cross covariances
% --------------------------------------------------------------------------
%%
%p.weighing=4; % MW weighing
%p.bias=1; % this deactivates an option I don't understand.

%p.weighing=1;
%p.bias=0;



% for testing purposes one can manually set these fields when p and
% branches are loaded from a dataset:
% fieldX = 'C6_mean_cycCor'; fieldY = 'muP9_fitNew_cycCor'; p.timeField ='C_time'
% fieldX = [FIELDPREFIX associatedFieldNames{2}]; fieldY=[FIELDPREFIX associatedFieldNames{3}]; p.timeField = [associatedFieldNames{1}]; p.bias=0; branches=branch_groups(1).branches; p.spacingError=0.05; p.weighing=1;

% Obtain the actual cross-correlations, variance and composite thereof:
% Note that the following settings will set the method to calculate the
% composite cross-correlation:
% p.weighing=1;     % Weighing method
% p.bias=0;         % Do not adjust 
% p.extraNorm=0;    % Normalization of separate branches
[branches_crossCov, crossCov_composite,total_weight] = MW_getCrossCov(p, branches, fieldX, fieldY, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm);
                  

%% determine composed cross-correlation


%%
% --------------------------------------------------------------------------
% Calc of correlations
% --------------------------------------------------------------------------
sourceField       = ['crossCov_' fieldX '_' fieldY];
sourceField_varX  = ['varX_' fieldX];
sourceField_varY  = ['varY_' fieldY];
targetField       = ['crossCor_' fieldX '_' fieldY];
%blubb % force truncation to 63 characters to suppress waring NW
%2012-05-16. AKWARD!!!!!
% edit MW 2017, let's try it without this truncation actually..
% if length(targetField)>63
%     targetField=targetField(1:63);
% end
% if length(sourceField_varX)>63
%     sourceField_varX=sourceField_varX(1:63);
%     warning('Field name too long, truncated..'); % MW line added
% end
% if length(sourceField_varY)>63
%     sourceField_varY=sourceField_varY(1:63);
%     warning('Field name too long, truncated..');
% end
% if length(sourceField)>63
%     sourceField=sourceField(1:63);
%     warning('Field name too long, truncated..');
% end
% blubb ende

% make branches
for br = 1:length(branches)
    
    % calculate normalization factor
    varXbr  = branches_crossCov(br).(sourceField_varX);
    varYbr  = branches_crossCov(br).(sourceField_varY);
    varXYbr = sqrt(varXbr*varYbr);
    
    %  normalize by variance 
    branches(br).(targetField) = branches_crossCov(br).(sourceField) / varXYbr;
    branches(br).(targetField) = branches_crossCov(br).(sourceField) / varXYbr;
    
    branches(br).corr_time = branches_crossCov(br).corr_time;
    branches(br).X = crossCov_composite.X;
    branches(br).Y = branches(br).(targetField);
end


% Normalize the composite covariance with variance to get the cross-correlation
% ===
crossCor_composite = struct;

% normalization factor
varX = crossCov_composite.(sourceField_varX); % variance in X 
varY = crossCov_composite.(sourceField_varY); % variance in Y
varXY = sqrt(varX*varY); % 

crossCor_composite.X = crossCov_composite.X;
crossCor_composite.Y = crossCov_composite.Y / varXY;

%% plotting new way

if isfield(p,'extraPlotting')
    %%
    figure(1); clf; hold on;

    % Plot separate branches' CCs
    % ===
    for idx=1:numel(branches)  
        plot(branches(idx).corr_time, branches(idx).(targetField),'-','LineWidth',1);
    end

    % Now plot the composite
    % ===
    plot(crossCor_composite.X, crossCor_composite.Y,'-o','LineWidth',3,'Color','k');

    % Cosmetics
    title(['R(' fieldX ',' 10 fieldY ')'],'Interpreter','None');
    ylabel('Correlation');
    xlabel('Lag');

    MW_makeplotlookbetter(20);

    myYlim=[-1,1];
    ylim(myYlim);
    plot([0,0],myYlim,'-k');

    myXlim=[min(branches(idx).corr_time), max(branches(idx).corr_time)];
    xlim(myXlim)
    plot(myXlim,[0,0],'-k');

end

% -------------------------------------------------------------------------
