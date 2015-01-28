% DJK_plotInFig_line plots line in current figure
% * not plotted when figu.LineWidth{} = 0
% * with data.show, for each point can be indicated whether it should be plotted (others replaced by NaN)
%
% REQUIRED INPUTS:
%  * data.X{}
%  * data.Y{}
%
% OPTIONAL INPUTS
%  * figu.LineStyle        - '-'
%  * figu.LineWidth        - 1
%  * figu.Color            - [0 0 0]
%  * figu.Marker           - 'none'
%  * data.show             - ones(size(data.X{1}))
%
% Code written by Daan Kiviet

function figu = DJK_plotInFig_line(data, figu);
% SET OPTIONAL INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' * ' mfilename]);
figu = DJK_setOptionalInput(figu, ['figu.LineStyle']       , '-');
figu = DJK_setOptionalInput(figu, ['figu.LineWidth']       , 1);
figu = DJK_setOptionalInput(figu, ['figu.Color']           , [0 0 0]);
figu = DJK_setOptionalInput(figu, ['figu.Marker']          , 'none');
data = DJK_setOptionalInput(data, ['data.show']            , ones(size(data.X{1})));

data = DJK_fillCellsToSize(data, length(data.X));
figu = DJK_fillCellsToSize(figu, length(data.X));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(data.X)
  if figu.LineWidth{i} > 0
    data.X{i}(find(~data.show{i})) = deal(NaN);
    data.Y{i}(find(~data.show{i})) = deal(NaN);
    line( data.X{i}, data.Y{i}, ...
          'LineStyle', figu.LineStyle{i}, ...
          'LineWidth', figu.LineWidth{i}, ...
          'Color', figu.Color{i}, ...
          'Marker', figu.Marker{i});
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
