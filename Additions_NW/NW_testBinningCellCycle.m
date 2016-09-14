function []=NW_testBinningCellCycle(p,schnitzcells,fieldY,numberBinArray,varargin)
% This function allows you to test several number of bins and interpolation
% methods, calculated via NW_plot_cellCycleDependence and compare them.
% -> Possibility to choose good numberBins for actual calculation
%
%
% OUTPUT
% plots
%
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzcells'  - schnitzcells structure, where useForPlot must be set
% 'fieldY'        - field from schnitzcells used for Y axis
% 'numberBinArray'  - Sets which bin numbers will be tested
%
%
% OPTIONAL ARGUMENTS:
% 'onScreen' = 0      automatically save and close image (default)
%              1      will only show image, and not ask to save
%                     NOTE: this refers to the mean and interpolation
%                     images. Single traces are always saved and closed
% 'DJK_saveDir'       Directory where images will be saved. Default:
%                     "p.analysisDir '\Cellcycle\'"
%                     Also Forwarded to NW_plotCellCycle
% 'selectionName'
% 'ylim'
% 'ylimCC'            Forwarded to NW_plotCellCycle
% 'interpolMethodCC'  Forwarded to NW_plotCellCycle
% 'periodicDataCC'    Forwarded to NW_plotCellCycle
% 'borderBinSizeCC'   Forwarded to NW_plotCellCycle
% 'equalBinsCC'       Forwarded to NW_plotCellCycle
% 'ignorePointsCC'    Forwarded to NW_plotCellCycle
% 'modeCC'            Forwarded to NW_plotCellCycle
% 'mytimefieldCC'     Forwarded to NW_plotCellCycle
%
% Note: more exotic features of NW_plotCellCycleDependence (time-reference-mode) 
% have to be changed directly here in this
% function and are not externally accessible

%akward
global extschnitzcounter;

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Copy 'p' to be able to parse it without extra fields to
% 'NW_plotCellCycleDependence'
p_intern=p;

% Settings
numRequiredArgs = 4; functionName = 'NW_testBinningCellCycle'; 

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
  p.onScreen = 0;
end
if ~existfield(p,'selectionName')
  p.selectionName = '';
end
if ~existfield(p,'interpolMethodCC') %special case
  p.interpolMethodCC = 'spline';
end

% Just in case it has not been set yet
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'cellCycle' filesep];
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
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Not plotting!');
  return;
end

if ~existfield(schnitzcells(1),fieldY)
  disp(['Field ' fieldY ' does not exist. Not plotting!']);
  return;
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Prepare for plotting
%--------------------------------------------------------------------------
%create colormap 
myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];

%x values for interpolation
xx=0:0.01:1;

% create text with input options for NW_plot_CellCycleDependence
myInput=['p_intern, schnitzcells, ''' fieldY ''', numbin, ''DJK_saveDir'', ''' p.DJK_saveDir ''' '];
if existfield(p,'ylimCC')
    myInput=[myInput, ', ''ylim'', [' num2str(p.ylimCC(1)) ' ' num2str(p.ylimCC(2)) '] ' ];
end
if existfield(p,'interpolMethodCC')
    myInput=[myInput, ', ''interpolMethod'', ''' p.interpolMethodCC ''''];
end
if existfield(p,'periodicDataCC')
    myInput=[myInput, ', ''periodicData'', ' num2str(p.periodicDataCC)];
end
if existfield(p,'borderBinSizeCC')
    myInput=[myInput, ', ''borderBinSize'', ' num2str(p.borderBinSizeCC)];
end
if existfield(p,'equalBinsCC')
    myInput=[myInput, ', ''equalBins'', ' num2str(p.equalBinsCC)];
end
if existfield(p,'ignorePointsCC')
    myInput=[myInput, ', ''ignorePoints'', ' num2str(p.ignorePointsCC)];
end
if existfield(p,'modeCC')
    myInput=[myInput, ', ''mode'', ''' p.modeCC ''''];
end
if existfield(p,'mytimefieldCC')
    myInput=[myInput, ', ''mytimefield'', ''' p.mytimefieldCC ''''];
end
%--------------------------------------------------------------------------
% Loop over numberBinArray and plot
%--------------------------------------------------------------------------
% figure names
fig1Name=['meanBinnedData_', fieldY];
fig2Name=['interpolation_', fieldY, '_', p.interpolMethodCC];

% initiate figures
fig1=figure(); clf;
scrsz = get(0, 'ScreenSize');
set(fig1,'Position', [50 50 scrsz(3)-200 scrsz(4)-150], 'Name', fig1Name);%, 'visible','off');
clf
hold on
fig2=figure(); clf;
scrsz = get(0, 'ScreenSize');
set(fig2,'Position', [50 50 scrsz(3)-200 scrsz(4)-150], 'Name', fig2Name);%, 'visible','off');
clf
hold on

for i=1:length(numberBinArray)
    numbin=numberBinArray(i);
    clear xdata ydata mypoly yy
    eval(['[mypoly,xdata,ydata]=NW_plot_cellCycleDependence(' myInput ');']);
    set(0,'CurrentFigure',fig1)
    plot(xdata,ydata,'.-','LineWidth',2,'Color',myColor(i,:),'MarkerSize',20)
       
    yy=ppval(mypoly,xx);
    set(0,'CurrentFigure',fig2)
    plot(xx,yy, '-','LineWidth',2,'Color',myColor(i,:))
    %really stupid solution to get sth to write a legend to
    if i==1
        for j=2:length(numberBinArray)
            plot(xdata(1),ydata(1),'LineWidth',2,'Color',myColor(j,:))
        end
        plot(xdata(1),ydata(1),'LineWidth',2,'Color',myColor(i,:))
    end     
    plot(xdata,ydata,'.','MarkerSize',15,'Color',myColor(i,:),'MarkerSize',20)
end

% create legend
mylegend=['''', num2str(numberBinArray(1)), '''' ];
for i=2:length(numberBinArray)
     mylegend=[mylegend, ', ''', num2str(numberBinArray(i)), ''''];
end
set(0,'CurrentFigure',fig1)
eval(['legend(', mylegend, ',''Location'',''S'')']);
title([p.movieDate ' ' p.movieName '.     ' fieldY '.   plotted ' num2str(extschnitzcounter) ' schnitzes'], ...
      'interpreter','none','FontWeight','bold','FontSize',12);
set(0,'CurrentFigure',fig2)  
eval(['legend(', mylegend, ',''Location'',''S'')']);
title([p.movieDate ' ' p.movieName '.     ' fieldY '.    ' p.interpolMethodCC, ' interpolation.   plotted ' num2str(extschnitzcounter) ' schnitzes'], ...
      'interpreter','none','FontWeight','bold','FontSize',12);

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% SAVE AND CLOSE FIGURE (if p.onScreen==0)
%--------------------------------------------------------------------------
if p.onScreen == 0
    saveSameSize(fig1,'file',[p.DJK_saveDir fig1Name '.png'], 'format', 'png');
    disp([' * Saved plot in ' fig1Name '.png']);
    saveSameSize(fig2,'file',[p.DJK_saveDir fig2Name '.png'], 'format', 'png');
    disp([' * Saved plot in ' fig2Name '.png']);
    close(fig1);
    close(fig2);
end
%--------------------------------------------------------------------------
