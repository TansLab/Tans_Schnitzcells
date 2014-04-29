
%% Settings

amolf_colors;

% Common settings
myRootDir = 'D:\MICROSCOPE_EXPERIMENTS\To_Analyze\'

%% Load data_______________________________________________________________

% 37 degrees data *********************************************************

p = DJK_initschnitz('pos1crop','2014-04-04','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

fitTime = [0 300];
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_37_pos1, ploty_37_pos1] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

% mutants
p = DJK_initschnitz('pos6crop','2014-04-04','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

fitTime = [0 300];
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_37_pos6, ploty_37_pos6] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);


% 42 degrees data *********************************************************

% Get data pos 2 and obtain scatter values_________________________________

p = DJK_initschnitz('pos2crop','2014-03-24','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

fitTime = [0 300];
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_pos2, ploty_pos2] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

% Get data pos 13 and obtain scatter values________________________________

p = DJK_initschnitz('pos13crop','2014-03-24','e.coli.AMOLF','rootDir',myRootDir, 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','none','fluor2','none','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

fitTime = [0 300];
s_all = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1); name_all = 'all';

[plotx_pos13, ploty_pos13] = ...
    DJK_plot_scatterColor(p, s_all, 'av_mu_fitNew', 'av_time', 'gen', 'ylim', [0 4], 'selectionName', name_all, 'plotRegression', 0, 'onScreen', 0);

%% Plotting the data________________________________________________________

% Scatter plot
% ====

close all;

figure(1);
hold on;

h2 = scatter(plotx_pos2, ploty_pos2, 100,colorAmolfGreen,'x','LineWidth',3); % WT, 42C
h13 = scatter(plotx_pos13, ploty_pos13, 100,colorAmolfBlue,'x','LineWidth',3); % 717, 42C
h42_pos1 = scatter(plotx_37_pos1, ploty_37_pos1, 100,colorAmolfGreen,'o','LineWidth',3); % WT, 37C
h42_pos6 = scatter(plotx_37_pos6, ploty_37_pos6, 100,colorAmolfBlue,'o','LineWidth',3); % 717, 37C


%set(h13)

set(gca,'FontSize',20);
title('Doubling time for each ''individual''');
xlabel('time (min)');
ylabel('\mu (doublings/hr)');

hold off;

%% make histograms
% =====

% bins
dx=2/10;
mybins=[0:dx:2];

% hist
[score2] = hist(ploty_pos2,mybins);
[score13] = hist(ploty_pos13,mybins);
[score_37_pos1] = hist(ploty_37_pos1,mybins);
[score_37_pos6] = hist(ploty_37_pos6,mybins);

% normalizing
A2=sum(score2)*dx;
A13=sum(score13)*dx;
A_37_pos1=sum(score_37_pos1)*dx;
A_37_pos6=sum(score_37_pos6)*dx;

score2=score2./A2;
score13=score13./A13;
score_37_pos1=score_37_pos1./A_37_pos1;
score_37_pos6=score_37_pos6./A_37_pos6;

%% plotting
figure(2);
hold on;

plot(mybins, score2,'-x','LineWidth',3,'color',colorAmolfGreen); % WT, 42C
plot(mybins, score13,'-x','LineWidth',3,'color',colorAmolfBlue); % 717, 42C
plot(mybins, score_37_pos1,'--o','LineWidth',3,'color',colorAmolfGreen); % WT, 37C
plot(mybins, score_37_pos6,'--o','LineWidth',3,'color',colorAmolfBlue); % WT, 37C

set(gca,'FontSize',20);
title('Distribution of doubling time');
xlabel('\mu (doublings/hr)');
ylabel('Probability (normalized)');


hold off;







