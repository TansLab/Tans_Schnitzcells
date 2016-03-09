% DJK_plotInFig_errorbar plots errorbars in current figure
% * not plotted when figu.LineWidth{} = 0
% * with data.show, for each point can be indicated whether it should be plotted (others replaced by NaN)
%
% REQUIRED INPUTS:
%  * data.X{}
%  * data.Y{}
%
% OPTIONAL INPUTS
%  * figu.LineStyle        - '-'
%  * figu.LineWidth        - 0.75
%  * figu.Color            - [0.5 0.5 0.5]
%  * figu.Marker           - 'none'
%  * figu.TeeLength        - 0.02 % lengths in units normalized relative to the x-axis.
%  * data.show             - ones(size(data.X{1}))
%  * data.X_error{}        - zeros(size(data.X{1}))
%  * data.Y_error{}        - zeros(size(data.X{1}))
%
% Code written by Daan Kiviet

function figu = DJK_plotInFig_errorbar(data, figu);
% SET OPTIONAL INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figu = DJK_setOptionalInput(figu, ['figu.LineStyle']       , '-');
figu = DJK_setOptionalInput(figu, ['figu.LineWidth']       , 0.75);
figu = DJK_setOptionalInput(figu, ['figu.Color']           , [0.5 0.5 0.5]);
figu = DJK_setOptionalInput(figu, ['figu.Marker']          , 'none');
figu = DJK_setOptionalInput(figu, ['figu.TeeLength']       , 0.02);
data = DJK_setOptionalInput(data, ['data.show']            , ones(size(data.X{1})));
data = DJK_setOptionalInput(data, ['data.X_error']         , zeros(size(data.X{1})));
data = DJK_setOptionalInput(data, ['data.Y_error']         , zeros(size(data.X{1})));

data = DJK_fillCellsToSize(data, length(data.X));
figu = DJK_fillCellsToSize(figu, length(data.X));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(data.X)
  if figu.LineWidth{i} > 0
    data.X{i}(find(~data.show{i}))        = deal(NaN);
    data.Y{i}(find(~data.show{i}))        = deal(NaN);
    data.X_error{i}(find(~data.show{i}))  = deal(NaN);
    data.Y_error{i}(find(~data.show{i}))  = deal(NaN);
    xtee                                  = figu.TeeLength{i} * (figu.axes.xlim(2)-figu.axes.xlim(1));
    ytee                                  = figu.TeeLength{i} * (figu.axes.ylim(2)-figu.axes.ylim(1)) * (figu.axes.width / figu.axes.height); 
    [errorBarX, errorBarY]                = DJK_getCoordinatesErrorbar(data.X{i}, data.Y{i}, data.X_error{i}, data.Y_error{i}, xtee, ytee);
    
    line( errorBarX, errorBarY, 'LineStyle', figu.LineStyle{i}, ...
                                'LineWidth', figu.LineWidth{i}, ...
                                'Color', figu.Color{i}, ...
                                'Marker', figu.Marker{i});
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errorBarX, errorBarY] = DJK_getCoordinatesErrorbar(dataX, dataY, errorX, errorY, xtee, ytee);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xl      = dataX - 0.5*xtee;
xr      = dataX + 0.5*xtee;
ytop    = dataY + errorY;
ybot    = dataY - errorY;

yt      = dataY + 0.5*ytee;
yb      = dataY - 0.5*ytee;
xleft   = dataX - errorX;
xright  = dataX + errorX;

% build up nan-separated vector for bars
errorBarX = NaN(length(dataX)*18,1);
errorBarY = NaN(length(dataY)*18,1);

if sum(errorY>0) > 0 % draw error in y
  % vertical line
  errorBarX( 1:18:end) = dataX;   errorBarY( 1:18:end) = ytop;
  errorBarX( 2:18:end) = dataX;   errorBarY( 2:18:end) = ybot;
  % top tee
  errorBarX( 4:18:end) = xl;      errorBarY( 4:18:end) = ytop;
  errorBarX( 5:18:end) = xr;      errorBarY( 5:18:end) = ytop;
  % bottom tee
  errorBarX( 7:18:end) = xl;      errorBarY( 7:18:end) = ybot;
  errorBarX( 8:18:end) = xr;      errorBarY( 8:18:end) = ybot;
end

if sum(errorX>0) > 0 % draw error in x
  % horizontal line
  errorBarX(10:18:end) = xleft;   errorBarY(10:18:end) = dataY;
  errorBarX(11:18:end) = xright;  errorBarY(11:18:end) = dataY;
  % left tee
  errorBarX(13:18:end) = xleft;   errorBarY(13:18:end) = yt;
  errorBarX(14:18:end) = xleft;   errorBarY(14:18:end) = yb;
  % right tee
  errorBarX(16:18:end) = xright;  errorBarY(16:18:end) = yt;
  errorBarX(17:18:end) = xright;  errorBarY(17:18:end) = yb;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%