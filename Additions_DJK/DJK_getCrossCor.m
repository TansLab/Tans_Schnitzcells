function [branches, crossCor_composite] = DJK_getCrossCor(p, branches, fieldX, fieldY, varargin);
% function [branches, crossCor_composite] = DJK_getCrossCor(p, branches, fieldX, fieldY, varargin);
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
%                 for larger delay times the autocorrelatiNature439_608on will always decrease
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
numRequiredArgs = 4; functionName = 'DJK_getCrossCor';

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
[branches_crossCov, crossCov_composite] = DJK_getCrossCov(p, branches, fieldX, fieldY, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm);

[branches_crossCov_XX, crossCov_XX_composite] = DJK_getCrossCov(p, branches, fieldX, fieldX, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm, ...
                                                          'rzeroonly', 1);

[branches_crossCov_YY, crossCov_YY_composite] = DJK_getCrossCov(p, branches, fieldY, fieldY, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm, ...
                                                          'rzeroonly', 1);
% --------------------------------------------------------------------------

%% MW - use to plot output if desired
%{
targetField='crossCov_noise_muP9_fitNew_cycCor_noise_muP9_fitNew_cycCor';
figure; clf; hold on;
%for idx=1:numel(branches)
%    plot(branches_crossCov_YY(idx).corr_time,branches_crossCov_YY(idx).(targetField),'-');
%end
halfwayindex=(numel(crossCov_YY_composite.Y)+1)/2;
sqrtvarxxvaryy=sqrt(crossCov_XX_composite.Y(halfwayindex)*crossCov_YY_composite.Y(halfwayindex));
plot(crossCov_composite.X,crossCov_composite.Y./sqrtvarxxvaryy,'k-','LineWidth',2)


myYlim=[-2,2];
plot([0,0],myYlim,'-k');

% checking
plot([-100,100],[-1,-1],'--k');
plot([-100,100],[1,1],'--k');

ylim(myYlim);

ylabel('Correlation');
xlabel('Lag');

MW_makeplotlookbetter(20);
%}

%%
% --------------------------------------------------------------------------
% Calc of correlations
% --------------------------------------------------------------------------
sourceField     = ['crossCov_' fieldX '_' fieldY];
sourceField_XX  = ['crossCov_' fieldX '_' fieldX];
sourceField_YY  = ['crossCov_' fieldY '_' fieldY];
targetField     = ['crossCor_' fieldX '_' fieldY];
%blubb % force truncation to 63 characters to suppress waring NW
%2012-05-16. AKWARD!!!!!
if length(targetField)>63
    targetField=targetField(1:63);
end
if length(sourceField_XX)>63
    sourceField_XX=sourceField_XX(1:63);
end
if length(sourceField_YY)>63
    sourceField_YY=sourceField_YY(1:63);
end
if length(sourceField)>63
    sourceField=sourceField(1:63);
end
% blubb ende

% make branches
for br = 1:length(branches)
  branches(br).(targetField) = branches_crossCov(br).(sourceField) / ...
        sqrt(branches_crossCov_XX(br).(sourceField_XX)(end) * branches_crossCov_YY(br).(sourceField_YY)(end));
  branches(br).corr_time = branches_crossCov(br).corr_time;
  branches(br).X = crossCov_composite.X;
  branches(br).Y = branches(br).(targetField);
end

% make composite
crossCor_composite = struct;
crossCor_composite.X = crossCov_composite.X;
r_0 = (length(crossCov_XX_composite.X) + 1) / 2; %=entry with time delay=0
crossCor_composite.Y = crossCov_composite.Y / sqrt(crossCov_XX_composite.Y(r_0)*crossCov_YY_composite.Y(r_0));

%% MW - use to plot output if desired
%{
figure; clf; hold on;
%for idx=1:numel(crossCor_composite)
plot(crossCor_composite(idx).X,crossCor_composite(idx).Y,'-');
%end

plot([min(crossCor_composite(idx).X),max(crossCor_composite(idx).X)],[0,0],'k-');
plot([0,0],[-1,1],'k-');

xlim([min(crossCor_composite(idx).X),max(crossCor_composite(idx).X)]);
ylim([-1,1]);
%}


% -------------------------------------------------------------------------
