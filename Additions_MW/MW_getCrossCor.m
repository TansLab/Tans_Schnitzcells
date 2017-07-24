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
%p.weighing=4; % MW weighing
%p.bias=1; % this deactivates an option I don't understand.

%p.weighing=1;
%p.bias=0;

p.weighing=1;
p.bias=0;
p.extraNorm=0;

[branches_crossCov, crossCov_composite,total_weight] = MW_getCrossCov(p, branches, fieldX, fieldY, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm);

% fieldX = 'noise_dC5_cycCor'; fieldY = 'muP9_fitNew_atdC5_cycCor'; p.timeField ='dC5_time'

% reset to before
p.weighing=3; 
p.bias=0; 

disp('Done');

%{
[branches_crossCov_XX, crossCov_XX_composite] = MW_getCrossCov(p, branches, fieldX, fieldX, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm, ...
                                                          'rzeroonly', 1);

[branches_crossCov_YY, crossCov_YY_composite] = MW_getCrossCov(p, branches, fieldY, fieldY, ...
                                                          'timeField', p.timeField, ...
                                                          'bias', p.bias, ...
                                                          'weighing', p.weighing, ...
                                                          'spacingError', p.spacingError, ...
                                                          'extraNorm', p.extraNorm, ...
                                                          'rzeroonly', 1);
%}

%
%figure; clf; hold on;
%plot(branches_crossCov(1).corr_time,branches_crossCov(1).w_crossCov_noise_C6_mean_cycCor_noise_muP9_fitNew_cycCor)
% --------------------------------------------------------------------------

%% plotting new way
targetField=['crossCov_' fieldX '_' fieldY];
%targetField='crossCov_noise_dC5_cycCor_noise_muP9_fitNew_atdC5_cycCor'
figure(1); clf; hold on;

% Plot separate branches' CCs
% ===
title(['R(' fieldX ',' fieldY ')'],'Interpreter','None');
halfwayidx=(numel(branches_crossCov(1).corr_time)+1)/2;
for idx=1:numel(branches)

    corrT{idx} = branches_crossCov(idx).corr_time;

    corrR{idx} = (branches_crossCov(idx).(targetField))./ ...
               ( ...
                sqrt(  branches_crossCov(idx).(['var_' fieldX]) * ...
                       branches_crossCov(idx).(['var_' fieldY])  ));
    
    %plot(corrT{idx}, corrR{idx},'-o','LineWidth',(numel(branches)-idx+1)*3);
    plot(corrT{idx}, corrR{idx},'-','LineWidth',1);
end

% Now calculate and plot the composite
% ===

varianceXXYY=sqrt(  crossCov_composite.(['var_' fieldX])./(total_weight(halfwayidx)) * ...
                    crossCov_composite.(['var_' fieldY])./(total_weight(halfwayidx))  );

plot(crossCov_composite.X, crossCov_composite.Y./varianceXXYY,'-o','LineWidth',3,'Color','k');

%plot(crossCov_composite.X, crossCov_composite.Y./varianceXXYY_noweight,'b--','LineWidth',2);



% checking
%plot([-100,100],[-1,-1],'--k');
%plot([-100,100],[1,1],'--k');


ylabel('Correlation');
xlabel('Lag');

MW_makeplotlookbetter(20);

myYlim=[-1,1];
ylim(myYlim);
plot([0,0],myYlim,'-k');

myXlim=[min(corrT{idx}), max(corrT{idx})];
xlim(myXlim)
plot(myXlim,[0,0],'-k');

%%
uniqueTimes=unique([corrT{1}]); 
warning('Times are not exactly identical -- check this'); % TODO!!!!
dT = uniqueTimes(2)-uniqueTimes(1);
edges = [uniqueTimes(1)-dT/2 uniqueTimes+dT/2];
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues, medianValuesForBins]=binnedaveraging(corrT,corrR,edges)
plot(binCenters,meanValuesForBins,'k-','LineWidth',3);

plot([min(edges),max(edges)],[0,0],'k-');
plot([0,0],[-1,1],'k-');
ylim([-1,1]);
xlim([min(edges),max(edges)]);

%%
total_weight
%var_weight = total_weight((length(total_weight)+1)/2)
%%

figure; clf; hold on;
title('calculating median for separately divided by variance');
corrTw={}; corrRw={};
for idx=1:numel(branches)

    corrTw{idx} = branches_crossCov(idx).corr_time;
    corrRw{idx} = branches_crossCov(idx).([targetField])./ ...
                   ( ...
                    sqrt(  branches_crossCov(idx).(['var_' fieldX]) * ...
                           branches_crossCov(idx).(['var_' fieldY])  ));

    plot(corrTw{idx},...
         corrRw{idx},...
         '-');
end

uniqueTimes=unique([corrTw{:}]); 
dT = uniqueTimes(2)-uniqueTimes(1);
edges = [uniqueTimes(1)-dT/2 uniqueTimes+dT/2];
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging(corrTw,corrRw,edges);
plot(binCenters,medianValuesForBins,'k-','LineWidth',3);


plot([min(edges),max(edges)],[0,0],'k-');
plot([0,0],[-1,1],'k-');
ylim([-1,1]);
xlim([min(edges),max(edges)]);



%% determine composed cross-correlation



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

%%
targetField='crossCov_noise_muP9_fitNew_cycCor_noise_muP9_fitNew_cycCor';
figure; clf; hold on;
for idx=1:numel(crossCor_composite)
    plot(crossCor_composite(idx).X,crossCor_composite(idx).Y,'-');
end

% -------------------------------------------------------------------------
