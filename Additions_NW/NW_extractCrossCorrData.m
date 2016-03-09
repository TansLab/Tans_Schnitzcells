function CorrData=DJK_plot_crosscorrelation_standard_error_store(p,branch_groups, fieldX, fieldY,varargin);
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
% Optional: CorrData: Matrix which contains: timepoint [min]. Crosscorr.
%                                            Normalized StdDev of CrossCorr
%                                            (StdDev/sqrt(#branchgroups)).
%                                            Crosscov
%                                            Normalized Stdev of Crosscov
%                     [t1, CC(t1), NoiseCC(t1), CrossCov(t1), NoiseCrossCov(t1);
%                     [t2, CC(t2), NoiseCC(t2), CrossCov(t2), NoiseCrossCov(t2);  ...]
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
% 'showPlots'     0: does not show any plots
%                 1: shows the plots (typically 'errornorm') which are
%                 specified in code (check below) -> default
% 
% %%% range of optional arguments could be extended %%%
%
% -------------------------------------------------------------------------

% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias = 0;      % 0 will adjust for less data at larger delay times
weighing = 2;  % 2 performs 3/4 weighing
extraNorm = 0; % 0 performs no extra normalization  %blubb

colorMode=1;   %=1 single branch groups are plotted in different colors. otherwise: plotted in gray
               % works only for up to 12 branches (by now)
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
if ~existfield(p,'showPlots ')
  p.showPlots = 1;
end


%--------------------------------------------------------------------------
% CALC COMPOSITE CORR FOR EACH GROUP 
% -------------------------------------------------------------------------
%p = struct; p.movieName = 'nonsense';  Noreen(11/2011): I commented it out,
%                                       otherwise the input-'p' to the whole function is overwritten, including
%                                       information about 'saveDir' (older version did not require
%                                       'p'as a function-input.)
for i = 1:length(branch_groups)
  [~, branch_groups(i).composite_corr, branch_groups(i).composite_cov ] = DJK_getCrossCorAndCov(p, branch_groups(i).branches, fieldX, fieldY, 'bias', bias, 'weighing', weighing, 'extraNorm', extraNorm);
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
  composite_cov(i,:) = branch_groups(i).composite_cov.Y;
  
end

% **** PLOT CROSSCORR ****
h=figure('visible','on'); 
% if single groups are to be plotted in color, get the colors
if colorMode==1 %start one color earlier than in NW_makemovieBranchgroups!
    myColor=[1 1 0.9; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
end
for i = 1:length(branch_groups)
    if colorMode==1
        colorindex=rem(i,size(myColor,1))+1;
          plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', myColor(colorindex,:)); hold on;
    else    
          plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
    end
end
plot(branch_groups(i).composite_corr.X/60, mean(composite_corr), '-', 'LineWidth', 3, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['single branch_groups (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
if colorMode==1 % legend with branch_group indices
    mylegend='''1'' ';
    for i=2:length(branch_groups)
        mylegend=[mylegend, ', ''', num2str(i),''''];
    end
    eval(['legend(' mylegend ')']);
end   
%saveas(gcf,[p.DJK_saveDir figureName1 '.png'], 'png');
close(h)

h=figure('visible','on'); 
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['errorbars (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
%saveas(gcf,[p.DJK_saveDir figureName2 '.png'], 'png');
close(h)

h=figure('visible','on'); 
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr)/sqrt(length(branch_groups)), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscorr');
title(['errorbars normalized (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
%saveas(gcf,[p.DJK_saveDir figureName3 '.png'], 'png');
if ~p.showPlots
    close(h)
end
% -------------------------------------------------------------------------

% *** PLOT CROSSCOV ***
h=figure('visible','on'); 
% if single groups are to be plotted in color, get the colors
if colorMode==1 %start one color earlier than in NW_makemovieBranchgroups!
    myColor=[1 1 0.9; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
end
for i = 1:length(branch_groups)
    if colorMode==1
        colorindex=rem(i,size(myColor,1))+1;
          plot(branch_groups(i).composite_cov.X/60, branch_groups(i).composite_cov.Y, '-', 'LineWidth', 2, 'Color', myColor(colorindex,:)); hold on;
    else    
          plot(branch_groups(i).composite_cov.X/60, branch_groups(i).composite_cov.Y, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
    end
end
plot(branch_groups(i).composite_cov.X/60, mean(composite_cov), '-', 'LineWidth', 3, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscov');
title(['single branch_groups (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
if colorMode==1 % legend with branch_group indices
    mylegend='''1'' ';
    for i=2:length(branch_groups)
        mylegend=[mylegend, ', ''', num2str(i),''''];
    end
    eval(['legend(' mylegend ')']);
end   
%saveas(gcf,[p.DJK_saveDir figureName1 '.png'], 'png');
close(h)

h=figure('visible','on'); 
errorbar(branch_groups(i).composite_cov.X/60,mean(composite_cov),std(composite_cov), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscov');
title(['errorbars (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
%saveas(gcf,[p.DJK_saveDir figureName2 '.png'], 'png');
close(h)

h=figure('visible','on'); 
errorbar(branch_groups(i).composite_cov.X/60,mean(composite_cov),std(composite_cov)/sqrt(length(branch_groups)), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
xlabel('time [h]');
ylabel('crosscov');
title(['errorbars normalized (' fieldX ', ' fieldY ')'],'interpreter','none');
hold on
x_lim=get(gca,'xlim');y_lim=get(gca,'ylim');
plot(x_lim,[0 0],'-k'); plot([0 0 ],y_lim,'-k');  % plot axis
%saveas(gcf,[p.DJK_saveDir figureName3 '.png'], 'png');
if ~p.showPlots
    close(h)
end
% -------------------------------------------------------------------------


clear CorrData;
CorrData(:,1)=branch_groups(i).composite_corr.X/60;
CorrData(:,2)=mean(composite_corr);
CorrData(:,3)=std(composite_corr)/sqrt(length(branch_groups));
CorrData(:,4)=mean(composite_cov);
CorrData(:,5)=std(composite_cov)/sqrt(length(branch_groups));



