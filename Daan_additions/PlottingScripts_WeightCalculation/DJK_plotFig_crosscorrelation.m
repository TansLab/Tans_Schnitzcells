% DJK_plotFig_crosscorrelation makes a figure with crosscorrelations
% * 1: plots X=0 and Y=0 lines
% * 2: plots weights
% * 3: plots errorbars
% * 4: plots crosscorrelations
%
% REQUIRED INPUTS:
%  * data.cc.X{}
%  * data.cc.Y{}
%
% OPTIONAL INPUTS:
%
% Code written by Daan Kiviet

function DJK_plotFig_crosscorrelation(data, figu);

% PLOT AND SAVE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figu = DJK_plotInFig_prepareFig(figu);

% plot X=0 line
if isfield(data, 'lineX0')
  DJK_plotInFig_line(data.lineX0, figu.lineX0);
end

% plot Y=0 line
if isfield(data, 'lineY0')
  DJK_plotInFig_line(data.lineY0, figu.lineY0);
end

% plot Weight lines
if isfield(data.cc, 'weight')
  DJK_plotInFig_line(data.cc.weight, figu.cc.weight);
end

% plot crosscorrelation errorbars
if isfield(data.cc, 'errorbar')
  DJK_plotInFig_errorbar(data.cc.errorbar, figu.cc.errorbar);
end

% plot crosscorrelations
DJK_plotInFig_line(data.cc, figu.cc);

DJK_plotInFig_showFig(figu);
DJK_plotInFig_saveFig(figu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
