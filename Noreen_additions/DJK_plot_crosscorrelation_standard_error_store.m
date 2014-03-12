function DJK_plot_crosscorrelation_standard_error_store(p,branch_groups, fieldX, fieldY,varargin);
% Daans function but saves the generated plots
%
% This function plots the cross-correlation of fieldX and fieldY. Further
% details of calculation can be found in DJK_getCrossCor and
% DJK_getCrossCov.
% 
% OUTPUT PLOTS:
% 1st Figure: Cross-correlations for each individual branch tree of the
%             branch_groups (i.e. all branches resulting from one initial cell.
%             The number of branch trees is given by the first input argument in
%             DJK_trim_branch_data, 'nr-branches'if the excel file is used
%             and is typically 4 or 8. 
%             Additionally the mean, so the complete cross-correlation
%             function, is plotted.
% 2nd Figure: Cross-correlation function with error-bars
% 3nd Figure: Cross-correlation function with error-bars. The size of the error-bars
%             here was divided by the square root of the number of
%             individual branch trees (the size of 'branch_groups')
%
%
% REQUIRED ARGUMENTS:
% 'p'
% 'branch_groups'
% 'fieldX'        field X of which crosscorrelation should be calculated
% 'fieldY'        field Y of which crosscorrelation should be calculated
%
%
% OPTIONAL ARGUMENTS:
% 'selectionName: Creates subfolder with this name in which plots will be
%                 stored
% 'DJK_saveDir'   directory where images will be saved. 
%                 default: [p.analysisDir 'branches\']
% 
% %%% range of optional arguments could be extended %%%
%
% -------------------------------------------------------------------------

% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias = 0;      % 0 will adjust for less data at larger delay times
weighing = 2;  % 2 performs 3/4 weighing
extraNorm = 0; % 0 performs no extra normalization  %blubb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'DJK_plot_crosscorrelation_standard_error_store';

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

% create SaveDir 
if ~existfield(p,'selectionName ')
  p.selectionName = ''
end
% Save directory
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'branches' filesep];
end
if ~existfield(p,'selectionName ')
  p.selectionName = '';
end
if length(p.selectionName) > 0
  p.DJK_saveDir = [p.DJK_saveDir p.selectionName filesep];
end  
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end



%--------------------------------------------------------------------------
% CALC COMPOSITE CORR FOR EACH GROUP 
% -------------------------------------------------------------------------
%p = struct; p.movieName = 'nonsense';  Noreen(11/2011): I commented it out,
%                                       otherwise the input-'p' to the whole function is overwritten, including
%                                       information about 'saveDir' (older version did not require
%                                       'p'as a function-input.)
for i = 1:length(branch_groups)
  [trash, branch_groups(i).composite_corr] = DJK_getCrossCor(p, branch_groups(i).branches, fieldX, fieldY, 'bias', bias, 'weighing', weighing, 'extraNorm', extraNorm);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% PLOTTING 
% -------------------------------------------------------------------------
% Make Figure Name
figureName1 = ['crosscorr_' fieldX ' _ ' fieldY 'single_branchgroups'];
figureName2 = ['crosscorr_' fieldX ' _ ' fieldY '_errors'];
figureName3 = ['crosscorr_' fieldX ' _ ' fieldY '_errors_norm'];


% actual plotting 
for i = 1:length(branch_groups)
  composite_corr(i,:) = branch_groups(i).composite_corr.Y;
end

figure; 
for i = 1:length(branch_groups)
  plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
end
plot(branch_groups(i).composite_corr.X/60, mean(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['single branch_groups (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
saveas(gcf,[p.DJK_saveDir figureName1 '.png'], 'png');

figure;
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['errorbars (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
saveas(gcf,[p.DJK_saveDir figureName2 '.png'], 'png');

figure;
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr)/sqrt(length(branch_groups)), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['errorbars normalized (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
saveas(gcf,[p.DJK_saveDir figureName3 '.png'], 'png');
% -------------------------------------------------------------------------


