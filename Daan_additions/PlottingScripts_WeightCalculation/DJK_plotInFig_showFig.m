% DJK_plotInFig_showFig shows figure depending on settings.
%
% REQUIRED INPUTS:
%  * figu.fig{}
%
% OPTIONAL INPUTS
%  * figu.figure.show        - 1
%
%
% Code written by Daan Kiviet

function DJK_plotInFig_showFig(figu);
% SET OPTIONAL INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' * ' mfilename]);
figu = DJK_setOptionalInput(figu, ['figu.figure.show']      , 1);

figu.figure = DJK_fillCellsToSize(figu.figure, length(figu.fig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(figu.fig)
  if figu.figure.show{i}
    set( figu.fig{i}, 'visible', 'on');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%