function DJK_plot_crosscorrelation_standard_error(branch_groups, fieldX, fieldY);

% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias = 0;      % 0 will adjust for less data at larger delay times
weighing = 2;  % 2 performs 3/4 weighing
extraNorm = 0; % 0 performs no extra normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALC COMPOSITE CORR FOR EACH GROUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = struct; p.movieName = 'nonsense';
for i = 1:length(branch_groups)
  [trash, branch_groups(i).composite_corr] = DJK_getCrossCor(p, branch_groups(i).branches, fieldX, fieldY, 'bias', bias, 'weighing', weighing, 'extraNorm', extraNorm);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(branch_groups)
  composite_corr(i,:) = branch_groups(i).composite_corr.Y;
end

figure; 
for i = 1:length(branch_groups)
  plot(branch_groups(i).composite_corr.X/60, branch_groups(i).composite_corr.Y, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
end
plot(branch_groups(i).composite_corr.X/60, mean(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;

figure;
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;

figure;
errorbar(branch_groups(i).composite_corr.X/60,mean(composite_corr),std(composite_corr)/sqrt(length(branch_groups)), '-', 'LineWidth', 2, 'Color', [0 0 0]); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
