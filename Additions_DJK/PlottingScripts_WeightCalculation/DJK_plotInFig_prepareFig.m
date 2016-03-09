% DJK_plotInFig_prepareFig prepares Figure for plotting, and returns figure settings with figure handle
%
% REQUIRED INPUTS:
%
% OPTIONAL INPUTS
%  * figu.figure.filename        - 'unnamed figure'
%  * figu.axes.xlim              - [-1 1]
%  * figu.axes.ylim              - [-1 1]
%  * figu.axes.width             - 5
%  * figu.axes.height            - 3
%  * figu.axes.lineWidth         - 1
%  * figu.axes.offset            - [1.0 1.0 0.2 0.2]
%  * figu.axes.yscale            - 'linear'
%  * figu.axes.xticks            - [-1:1]
%  * figu.axes.yticks            - [-1:1]
%  * figu.axes.ticklength        - 0.015
%  * figu.axes.tick_fontsize     - 7
%  * figu.axes.label_fontsize    - 8
%  * figu.axes.xlabel            - 'unknown'
%  * figu.axes.ylabel            - 'unknown'
%  * figu.axes.xlabel_offset     - [0.0 0.0]
%  * figu.axes.ylabel_offset     - [0.0 0.0]
%
% Code written by Daan Kiviet

function figu = DJK_plotInFig_prepareFig(figu);
% SET OPTIONAL INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' * ' mfilename]);
figu = DJK_setOptionalInput(figu, ['figu.figure.filename']      , 'unnamed figure');
figu = DJK_setOptionalInput(figu, ['figu.axes.xlim']            , [-1 1]);
figu = DJK_setOptionalInput(figu, ['figu.axes.ylim']            , [-1 1]);
figu = DJK_setOptionalInput(figu, ['figu.axes.width']           , 5);
figu = DJK_setOptionalInput(figu, ['figu.axes.height']          , 3);
figu = DJK_setOptionalInput(figu, ['figu.axes.lineWidth']       , 1);
figu = DJK_setOptionalInput(figu, ['figu.axes.offset']          , [1.0 1.0 0.2 0.2]);
figu = DJK_setOptionalInput(figu, ['figu.axes.yscale']          , 'linear');
figu = DJK_setOptionalInput(figu, ['figu.axes.xticks']          , [-1:1]);
figu = DJK_setOptionalInput(figu, ['figu.axes.yticks']          , [-1:1]);
figu = DJK_setOptionalInput(figu, ['figu.axes.ticklength']      , 0.015);
figu = DJK_setOptionalInput(figu, ['figu.axes.tick_fontsize']   , 7);
figu = DJK_setOptionalInput(figu, ['figu.axes.label_fontsize']  , 8);
figu = DJK_setOptionalInput(figu, ['figu.axes.xlabel']          , 'unknown');
figu = DJK_setOptionalInput(figu, ['figu.axes.ylabel']          , 'unknown');
figu = DJK_setOptionalInput(figu, ['figu.axes.xlabel_offset']   , [0.0 0.0]);
figu = DJK_setOptionalInput(figu, ['figu.axes.ylabel_offset']   , [0.0 0.0]);

figu.figure = DJK_makeFieldsCells(figu.figure);
figu.axes = DJK_makeFieldsCells(figu.axes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PREPARE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0, 'ScreenSize'); scrwidth = scrsz(3); scrheight = scrsz(4); % make figure with close to full screen size, but hidden
dFBorder = 100; % distance From Border

for i = 1:length(figu.figure.filename)
  % CREATE FIGURE
  figu.fig{i} = figure( 'Position', [dFBorder dFBorder scrwidth-2*dFBorder scrheight-2*dFBorder], ...
                        'Name', figu.figure.filename{i}, 'visible','off'); % Position / screenSize: [left bottom width height]  

  % SET BASIC FIGURE SETTINGS
  set(figu.fig{i}, 'Units', 'centimeter');
  title('');
  legend('off'); % first remove legend, as it gives an extra axis
  axis([figu.axes.xlim{i} figu.axes.ylim{i}]);

  % GET AXIS HANDLES
  figAx            = get(figu.fig{i},'Children');
  figAxXLabel      = get(figAx,'XLabel');
  figAxYLabel      = get(figAx,'YLabel');

  % CHANGE AXES -> SIZE
  set(figAx, 'Units',     'centimeter');
  set(figAx, 'Box',       'on');
  set(figAx, 'LineWidth', figu.axes.lineWidth{i});
  set(figAx, 'Position',  [figu.axes.offset{i}(2) figu.axes.offset{i}(1) figu.axes.width{i} figu.axes.height{i}]);
  set(figAx, 'XLim',      figu.axes.xlim{i}, 'YLim', figu.axes.ylim{i});
  set(figAx, 'YScale',    figu.axes.yscale{i});

  % CHANGE AXES -> TICKS
  set(figAx,'XTick', figu.axes.xticks{i}, 'YTick', figu.axes.yticks{i}, ...
            'ticklength', [figu.axes.ticklength{i} figu.axes.ticklength{i}], 'XMinorTick', 'off', 'YMinorTick', 'off'); 

  % CHANGE AXES -> LABELS
  set(figAx,       'FontName', 'Arial', 'FontWeight', 'normal', 'FontSize', figu.axes.tick_fontsize{i});
  set(figAxXLabel, 'FontName', 'Arial', 'FontWeight', 'normal', 'FontSize', figu.axes.label_fontsize{i}, ...
                   'String', figu.axes.xlabel{i}, 'Units', 'centimeter');
  set(figAxYLabel, 'FontName', 'Arial', 'FontWeight', 'normal', 'FontSize', figu.axes.label_fontsize{i}, ...
                   'String', figu.axes.ylabel{i}, 'Units', 'centimeter');
  set(figAxXLabel, 'Position', [figu.axes.xlabel_offset{i} 0] + get(figAxXLabel,'Position'));
  set(figAxYLabel, 'Position', [figu.axes.ylabel_offset{i} 0] + get(figAxYLabel,'Position'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%