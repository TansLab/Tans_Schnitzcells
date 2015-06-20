% DJK_plot_scatterColor makes a scatter plot of particular data in
% schnitzcells. If schnitzcells has the field 'useForPlot', only schnitzes
% that are 1 there, will be plotted, otherwise all cells with data will be
% plotted. If fieldColor='none', no color will be used, and no colorbar
%
% OUTPUT
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzcells'  - schnitzcells structure, where useForPlot must be set
% 'fieldY'        - field from schnitzcells used for Y axis
% 'fieldX'        - field from schnitzcells used for X axis
% 'fieldColor'    - field from schnitzcells used for color ('none' if not used)
%
% OPTIONAL ARGUMENTS:
% 'colorNormalize'    Normalization values for color. If not given, will
%                     use min and max
%
% 'onScreen' = 0      automatically save and close image
%              1      will ask to save and close image (default)
%              2      will only show image, and not ask to save
% 'DJK_saveDir'       Directory where images will be saved. Default:
%                     "p.analysisDir 'schnitzcells\'"
% 'markerSize'        default: 8, used to be 20
% 'plotRegression'    default=1, else does not show
% 'selectionName'
% 'xlim'
% 'ylim'
% 'fieldColor_func'   default: first value -> @(x) x(1)
%                     alternative: min -> @(x) min(x)
%                     alternative: max -> @(x) max(x)
%                     alternative: mean -> @(x) mean(x)
% 'useArrayData'
%

function [plot_X, plot_Y] = DJK_plot_scatterColor(p, schnitzcells, fieldY, fieldX, fieldColor, varargin);


%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 5; functionName = 'DJK_plot_scatterColor';

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
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 1;
end
% This normalization settings are used
if ~existfield(p,'colorNormalize')
  colorNormalize = false; % will set based on fieldColor values later
else
  colorNormalize = true; 
end
if ~existfield(p,'fieldColor_func')
  p.fieldColor_func = @(x) x(1);
end
if ~existfield(p,'markerSize ')
  p.markerSize = 8;
end
if ~existfield(p,'plotRegression ')
  p.plotRegression = 1;
end
if ~existfield(p,'selectionName ')
  p.selectionName = '';
end
% Just in case it has not been set yet
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzcells' filesep];
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



%--------------------------------------------------------------------------
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Not plotting!');
  return;
end

if ~existfield(schnitzcells(1),fieldX)
  disp(['Field ' fieldX ' does not exist. Not plotting!']);
  return;
end

if ~existfield(schnitzcells(1),fieldY)
  disp(['Field ' fieldY ' does not exist. Not plotting!']);
  return;
end

useColor = 1;
if ~existfield(schnitzcells(1),fieldColor)
  if strcmp('NONE',upper(fieldColor))
    useColor = 0;
  else
    disp(['Field ' fieldColor ' does not exist. Not plotting!']);
    return;
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get data to plot
%--------------------------------------------------------------------------
useAllcells = ~existfield(schnitzcells(1),'useForPlot');
  
count = 0;
for cell = 1:length(schnitzcells)
  if useAllcells | schnitzcells(cell).useForPlot
    if length(schnitzcells(cell).(fieldX)) & length(schnitzcells(cell).(fieldY)) & (~useColor | length(schnitzcells(cell).(fieldColor)))
      if ~isnan(schnitzcells(cell).(fieldX)(1)) & ~isnan(schnitzcells(cell).(fieldY)(1)) % remove NaN values

          % Add single values to variables to plot plot_X, plot_Y (plot_Z
          % encodes color) -  idx_schnitzarray loops over 
          % vectors that are contained in the Schnitz.
          for idx_schnitzarray = 1:length(schnitzcells(cell).(fieldX)) % 2014/6/5 MW
              
                count = count + 1;
                plot_X(count) = schnitzcells(cell).(fieldX)(idx_schnitzarray);
                plot_Y(count) = schnitzcells(cell).(fieldY)(idx_schnitzarray);
                if useColor
                  plot_Z(count) = p.fieldColor_func(schnitzcells(cell).(fieldColor));
                end 
                
          end
          
      end
    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Regression analysis on data
%--------------------------------------------------------------------------
[stat_R,stat_P] = corrcoef(plot_X,plot_Y);
stat_R = stat_R(2,1);
stat_P = stat_P(2,1);
stat_fitCoef1 = polyfit(plot_X,plot_Y,1);

stat_stdX = std(plot_X);
stat_stdY = std(plot_Y);
stat_meanX = mean(plot_X);
stat_meanY = mean(plot_Y);
stat_n = length(plot_X);
stat_syx = sqrt( (stat_n-1) / (stat_n-2) * (stat_stdY*stat_stdY - stat_fitCoef1(1)*stat_fitCoef1(1)*stat_stdX*stat_stdX));
stat_sb = 1 / sqrt(stat_n-1) * stat_syx / stat_stdX;
stat_ta = DJKcopy_tinv(1-0.05/2,stat_n-2);
stat_plusmin = stat_ta * stat_sb;
stat_rxy_2 = 1 - (stat_syx*stat_syx)/(stat_stdY*stat_stdY);
% disp(['stat_R   : ' num2str(stat_R)]);
% disp(['stat_P   : ' num2str(stat_P)]);
% disp(['fit line : y = ' num2str(stat_fitCoef1(1)) ' * x + ' num2str(stat_fitCoef1(2))]);
% disp([num2str(stat_rxy_2) ' of variance in Y is accounted for by X']); 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Scale Color data
%--------------------------------------------------------------------------
if useColor
  number_of_colors = 100;
  map = colorGray(number_of_colors);

  if ~colorNormalize % was not set
    p.colorNormalize = [min(plot_Z) max(plot_Z)];
  end

  plot_Z_scaled = round( DJK_scaleRange(plot_Z, [p.colorNormalize(1) p.colorNormalize(2)], [1 number_of_colors]) );
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------
% Make Figure Name
figureName = [fieldY ' ___ ' fieldX];
figureFileName = [fieldY '___' fieldX];
if useColor
  figureName = [figureName ' ___ color: ' fieldColor];
  figureFileName = [figureFileName '___' fieldColor];
end

% make figure with full screen size
scrsz = get(0, 'ScreenSize');
% fig1 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)], 'Name', figureName);
fig1 = figure('Position', [100 100 scrsz(3)-125 scrsz(4)-200], 'Name', figureName, 'visible','off');
hold on;

if useColor
  % loop over colors, and plot points with that color
  for k=1:size(map,1) 
  %for k=1:size(map,1)*0.3
    if any(plot_Z_scaled==k)
      line( 'Xdata',plot_X(plot_Z_scaled==k), ...
            'Ydata',plot_Y(plot_Z_scaled==k), ...
            'LineStyle','none', ...
            'Color',map(k,:), ...
            'Marker','.', ...
            'MarkerSize',p.markerSize);
    end
  end
else
  line( 'Xdata',plot_X, ...
        'Ydata',plot_Y, ...
        'LineStyle','none', ...
        'Color','k', ...
        'Marker','.', ...
        'MarkerSize',p.markerSize);
end  
  
if p.plotRegression
  % add fitted line
  hold on;
  stat_fittedY = stat_fitCoef1(1)*plot_X + stat_fitCoef1(2);
  line( 'Xdata',plot_X, ...
        'Ydata',stat_fittedY, ...
        'LineStyle','-', ...
        'LineWidth',2, ...
        'Color','k', ...
        'Marker','none');

  % add statistics text
  text(0.02,0.98,['Correlation Coef R / R2 : ' DJK_setDecimalPlaces(stat_R,2) ' / ' DJK_setDecimalPlaces(stat_R*stat_R,2)],      'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
  text(0.02,0.96,['Prob R could be 0 (P)   : ' DJK_setDecimalPlaces(stat_P,4)],                                   'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
  text(0.02,0.93,['Regression line         : y = ' DJK_setDecimalPlaces(stat_fitCoef1(1),3) ' * x + ' DJK_setDecimalPlaces(stat_fitCoef1(2),2)], 'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
  text(0.02,0.91,['P=0.05 slope is between : ' DJK_setDecimalPlaces(stat_fitCoef1(1) - stat_plusmin,3) ' & ' DJK_setDecimalPlaces(stat_fitCoef1(1) + stat_plusmin,3)], 'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);

  text(0.02,0.89,['(Assuming X contains no observational errors.)'],      'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11); % added MW 2014/06/05
  disp('WARNING: Least squares assumes X contains no observational errors (noise).'); % added MW 2014/06/05
  disp('   Use total least squares (TLS) if X and Y both contain observational errors.'); % added MW 2014/06/05
  disp('   (This functionality is not available in this code.)'); % added MW 2014/06/05
  
  text(0.02,0.03,['var in Y caused by X    : ' DJK_setDecimalPlaces(stat_rxy_2,2)], 'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);

  text(0.80,0.05,['noise Y : ' DJK_setDecimalPlaces(stat_stdY/stat_meanY,2)], 'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
  text(0.80,0.03,['noise X : ' DJK_setDecimalPlaces(stat_stdX/stat_meanX,2)], 'sc','FontName','FixedWidth','FontWeight','bold','FontSize',11);
  
end

% label axes
xlabel(fieldX,'interpreter','none','FontWeight','bold','FontSize',12); % interpreter to avoid problem with underscores
ylabel(fieldY,'interpreter','none','FontWeight','bold','FontSize',12);

% in case xlim or ylim
if existfield(p,'xlim ')
  xlim(p.xlim);
end
if existfield(p,'ylim ')
  ylim(p.ylim);
end

% Add title
title([p.movieDate ' ' p.movieName ' for ' num2str(length(plot_X)) ' schnitzes in folder: schnitzcells\\' p.selectionName], ...
      'interpreter','none','FontWeight','bold','FontSize',12);

% % resize
% set(gca,'Position', [0.1 0.1 0.8 0.8]); %[0.13 0.11 0.69422 0.815]
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% ADD COLOR BAR
%--------------------------------------------------------------------------
if useColor
  colormap(map);
  cBar = colorbar('EastOutside'); %, 'peer',axes_handle

  % lowest value
  labelText = {p.colorNormalize(1)}; 
  labelHeight = [1];

  % min color value
  minColorHeight = round( DJK_scaleRange(min(plot_Z), [p.colorNormalize(1) p.colorNormalize(2)], [1 1+number_of_colors]) );
  if minColorHeight>1 & minColorHeight<1+number_of_colors
    labelText = [labelText 'min']; 
    labelHeight = [labelHeight minColorHeight]; 
  end  

  % max color value
  maxColorHeight = round( DJK_scaleRange(max(plot_Z), [p.colorNormalize(1) p.colorNormalize(2)], [1 1+number_of_colors]) );
  if maxColorHeight>1 & maxColorHeight>minColorHeight & maxColorHeight<1+number_of_colors
    labelText = [labelText 'max']; 
    labelHeight = [labelHeight maxColorHeight]; 
  end  

  % highest value
  labelText = [labelText p.colorNormalize(2)]; 
  labelHeight = [labelHeight 1+number_of_colors]; 

  % set Ticks
  set(cBar, 'YTickLabel', labelText, 'YTick', labelHeight);

  % set Label at color bar
  set(get(cBar,'YLabel'),'String',fieldColor,'interpreter','none','FontWeight','bold','FontSize',12);

  % reduce width colorbar
  graphPos = get(gca,'Position');
  cBarPos = get(cBar,'Position');
  cBarPos(3) = 0.5*cBarPos(3);
  set(cBar,'Position', cBarPos);
  set(gca,'Position', graphPos); % must do this, else wrong
end
%--------------------------------------------------------------------------

        
%--------------------------------------------------------------------------
% ASKING TO SAVE FIGURE
%--------------------------------------------------------------------------
if p.onScreen == 2
  set(fig1,'visible','on');
else
  % Ask to save the figure
  if p.onScreen
    set(fig1,'visible','on');
    saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
    pause(0.2);
  else
    saveFigInput='Yes and Close';
  end
  if (upper(saveFigInput(1))=='Y')
  %     saveas(fig1,[p.DJK_saveDir figureFileName '.fig']);
      saveSameSize(fig1,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');
      if (strcmp(saveFigInput,'Yes and Close'))
          close(fig1);
          pause(0.2);
      end
    disp([' * Saved plot in ' p.DJK_saveDir figureFileName '.png']);
  end
end
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function number_string = DJK_setDecimalPlaces(number, decimal_places);
number_string = sprintf(['%1.' num2str(decimal_places) 'f'], number);

% factor = power(10,decimal_places);
% number = number*factor;
% number = floor(number);
% number = number / factor;
% number_string = num2str(number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
