% DJK_makeFig_crosscorrelationWeight calculates weight of crosscorrelation
% and makes a figure of this part of crosscorrelation of E with mu and p with mu 
%
% Code written by Daan Kiviet

function DJK_makeFig_crosscorrelationWeight(varargin);

% DATA SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = struct;
dt.cc.expName                 = {'ICD' 'ICD' };
dt.cc.fitTime                 = {[102 2198] [102 2198]};
dt.cc.showRange               = {[-6 6] [-6 6]};
dt.cc.dataField1              = {'Y6_mean' 'dY5PN'};
dt.cc.dataField2              = {'muP23_fitNew' 'muP23_fitNew'};
dt.cc.timeField               = {'Y_time' 'dY5PN_time'};
dt.cc.cycleCorrect            = 1;
dt.cc.bias                    = 0;
dt.cc.weighing                = 2;
dt.cc.extraNorm               = 0;
dt.cc.groups                  = 4;

dt.cc.errorbar.barType        = {'std' 'std'};
dt.cc.errorbar.showSelection  = {[1 0] [1 0]};

dt.cc.errorbar.Pcutoff        = 1.00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SAVE/SHOW SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg = struct;
fg.figure.show                = 1; % 0 (do not show), 1 (show)
fg.figure.save                = 1; % 0 (no save), 1 (save), 2 (ask to save)
fg.figure.filetypes           = { {'pdf'} };
fg.figure.saveInExpDir        = 0;
fg.figure.foldername          = [''];
fg.figure.filename            = ['crosscorrelationWeight'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DATA PLOTTING PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg.cc.LineWidth               = {2 1}; 
fg.cc.Color                   = {[159 135 191]/256 [116  78 160]/256};

fg.cc.errorbar.LineWidth      = 0.75;
fg.cc.errorbar.Color          = DJK_lightenColors(fg.cc.Color);
fg.cc.errorbar.teeLength      = 0.02; % lengths in units normalized relative to the x-axis.

fg.cc.weight.LineWidth        = fg.cc.LineWidth; % 0 does not show weight
fg.cc.weight.Color            = fg.cc.Color;
fg.cc.weight.yoffset          = {0.05 0.10};

fg.lineX0.LineWidth           = 1;   % 0 will not draw it
fg.lineX0.Color               = [0 0 0];
fg.lineY0.LineWidth           = 1;   % 0 will not draw it
fg.lineY0.Color               = [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FIGURE AND AXES SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg.axes.width                 = 3.55; % width of axes in cm
fg.axes.height                = 2.13; % height of axes in cm
fg.axes.offset                = [0.9 1.0 0.2 0.2]; % offset between axis and figure in cm [bottom left top right]

fg.axes.xlabel                = 'time delay (h)';
fg.axes.ylabel                = 'cross-correlation';

fg.axes.xlabel_offset         = [ 0.0  0.1];
fg.axes.ylabel_offset         = [ 0.2  0.0];

fg.axes.xlim                  = [-6 6];
fg.axes.ylim                  = [-0.1 0.2];
fg.axes.xticks                = [-6:2:6];
fg.axes.yticks                = [-0.1:0.1:0.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = load(['ICD_RESULT_merged_crosscorrs.mat']);
dt.cc.errorbar.timedelay{1}        = data.time_merge_conc';
dt.cc.errorbar.timedelay{2}        = data.time_merge_rate';
dt.cc.errorbar.correlation{1}      = data.composite_corr_merge_conc_mu;
dt.cc.errorbar.correlation{2}      = data.composite_corr_merge_rate_mu;
for i=1:2
    dt.cc.errorbar.meanCorrelation{i}  = mean(dt.cc.errorbar.correlation{i});
    dt.cc.errorbar.STD{i}              = std(dt.cc.errorbar.correlation{i});
    [hypothesis{i}, dt.cc.errorbar.statP{i}] = ttest( dt.cc.errorbar.correlation{i} );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% PREPARE ERROR BARS AND WEIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(dt.cc.expName)
  
  % SELECT TO SHOW ERROR BAR TYPE
  switch dt.cc.errorbar.barType{i}
    case 'no'
      dt.cc.X{i}                = dt.cc.timedelay{i};
      dt.cc.Y{i}                = dt.cc.correlation{i};
      dt.cc.errorbar.Y_error{i} = zeros(size(dt.cc.Y{i}));
    case 'std'
      dt.cc.X{i}                = dt.cc.errorbar.timedelay{i};
      dt.cc.Y{i}                = dt.cc.errorbar.meanCorrelation{i};
      dt.cc.errorbar.Y_error{i} = dt.cc.errorbar.STD{i};
    case 'ste'
      dt.cc.X{i}                = dt.cc.errorbar.timedelay{i};
      dt.cc.Y{i}                = dt.cc.errorbar.meanCorrelation{i};
      dt.cc.errorbar.Y_error{i} = dt.cc.errorbar.STE{i};
  end
  dt.cc.errorbar.X{i}       = dt.cc.X{i};
  dt.cc.errorbar.Y{i}       = dt.cc.Y{i};
  dt.cc.errorbar.X_error{i} = zeros(size(dt.cc.errorbar.Y_error{i}));

  % SET SHOWRANGE OF CORRELATION
  dt.cc.show{i}               = zeros(size(dt.cc.X{i})); 
  idx                         = find(dt.cc.X{i}>=dt.cc.showRange{i}(1) & dt.cc.X{i}<=dt.cc.showRange{i}(2));
  dt.cc.show{i}(idx)          = deal(1); 

  % SET SHOWRANGE AND SHOWSELECTION OF ERRORBARS
  dt.cc.errorbar.show{i}      = zeros(size(dt.cc.errorbar.X{i})); 
  idx                         = [rem(floor(length(dt.cc.errorbar.X{i})/2)+dt.cc.errorbar.showSelection{i}(2),dt.cc.errorbar.showSelection{i}(1))+1 : dt.cc.errorbar.showSelection{i}(1) : length(dt.cc.errorbar.X{i})];
  dt.cc.errorbar.show{i}(idx) = deal(1); 
  idx                         = find(dt.cc.errorbar.X{i}<dt.cc.showRange{i}(1) | dt.cc.errorbar.X{i}>dt.cc.showRange{i}(2));
  dt.cc.errorbar.show{i}(idx) = deal(0); 
  
  % SET RANGE FOR SIGNIFICANCE
  idx                         = find(dt.cc.errorbar.statP{i}>=dt.cc.errorbar.Pcutoff);
  dt.cc.Y{i}(idx)             = deal(0);
  dt.cc.errorbar.show{i}(idx) = deal(0);
  
  % CALC WEIGHT FOR PLOTTED DATA
  dataX                       = dt.cc.X{i}(find(dt.cc.show{i}));
  dataY                       = dt.cc.Y{i}(find(dt.cc.show{i}));
  dt.cc.corWeight{i}          = sum( (dataX .* dataY) / sum(dataY) );
  dt.cc.weight.X{i}           = [dt.cc.corWeight{i} 0];
  dt.cc.weight.Y{i}           = [fg.axes.ylim(2)-fg.cc.weight.yoffset{i}*(fg.axes.ylim(2)-fg.axes.ylim(1)) fg.axes.ylim(2)-fg.cc.weight.yoffset{i}*(fg.axes.ylim(2)-fg.axes.ylim(1))];
  disp([num2str(i) ': Weight of ' dt.cc.expName{i} ' = ' num2str(dt.cc.corWeight{i})]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% GET Y=0 & X=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt.lineX0.X{1} = [0 0];
dt.lineX0.Y{1} = fg.axes.ylim;
dt.lineY0.X{1} = fg.axes.xlim;
dt.lineY0.Y{1} = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fg.cc.errorbar.axes     = fg.axes; % needed to determine tee length


DJK_plotFig_crosscorrelation(dt, fg);