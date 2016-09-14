function []=NW_plot_dependence_on_position(p,plotField,varargin)
% Plots the spatial dependence of plotField on the position of the schnitz
% within the colony. The function was originally written to investigate
% whether central cells within a colony grow more slowly (nutrient
% depletion etc.), however it also works for a fluorescence color as
% input in 'plotField'.
%
% Two plots are generated:
% 1) For each frame an image similar to a segementation image is made.
%    Cells are colored according to their value of 'plotField'.  
% 2) Correlation plot between 'plotField' and position. For position, the
%    distance to the centre of the colony is taken.
%
% Plots can be turned off within the function by changing show2dimplot,
% showscatter_plot_arithm/showscatter_plot_mass.
%
% Note: Every schnitz in one frame is equally weighted, independent of which
% generation it belongs to.
%
% REQUIRED ARGUMENTS:
% 'p'
% 'plotField'      field for which dependence on position in colony is
%                   plotted. Will typically be muPxx_fitNew_all (mu), or
%                   Y6_mean_all, fitted_Y6_mean, noise_Y6_mean etc
%                   (fluorescence colors) (if value=Nan for certain frame,
%                   image will stay empty)
%
% OPTIONAL ARGUMENTS:
% 'frameRange'      Frame numbers for which plots are calculated. Per
%                   default: framenr=300 
% 'colorRange'      Sets borders of colorbar (what input value will be
%                   plotted in red resp. green).
%                   If input is 'automatic', range will be calculated
%                   automatically from the mean and stddev of plotField of
%                   ALL schnitzes (including sick ones, excluded ones,
%                   first frames,...)
%                   Range will be set to [mean-colorWidth*stddev,
%                   mean+colorWidth*stddev]. colorWidth=2 is a good
%                   prefactor but can be changed within the program.
%
% 'myColorMap'      Which colormap to choose. By default will be 'hot'. 
%                   Usually a red-green colormap will be loaded before
%                   calling this function.
% NW_saveDir        saveDirectory. per default:
%                   analysisDir/schnitzells/spatialDependence
% rm_SchnitzNrs      remove Schnitzes from plotting (because they are also
%                   removed in analysis). default=[]


%--------------------------------------------------------------------------
% For debugging. turn plotting on/off
%--------------------------------------------------------------------------
show2dimplot=1;           % (1) segmentation-like image of cell colony
showscatterplot_arithm=1; % (2) colony centre: arithmetic mean of cell positions (maybe more mathematically correct)
showscatterplot_mass=0;   % (2) colony centre: centre of mass (greater influence of longer cells. maybe more biological)
                          % differences are only small after colony size of
                          % ca. 50 (very sufficient!) is reached
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'NW_plot_dependence_on_position';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
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

% Get frame numbers which will be plotted
if ~existfield(p,'frameRange')
  p.frameRange=300; 
end
% set ColorRange
if ~existfield(p,'colorRange')
  p.colorRange=[0.3 0.7];
end
% Just in case it has not been set yet
if ~existfield(p,'NW_saveDir')
  p.NW_saveDir = [p.analysisDir  'spatialDependence' filesep];
end
% Set colormap for first plot
if ~existfield(p,'myColorMap')
  p.myColorMap = colormap('hot');
  p.myColorMap(1,:)=[0.7 0.7 0.7];
end
% Make sure that DJK_saveDir directory exists
if exist(p.NW_saveDir)~=7
  [status,msg,id] = mymkdir([p.NW_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.NW_saveDir ' : ' msg]);
    return;
  end
end
% Set rm_SchnitzNrs
if ~existfield(p,'rm_SchnitzNrs')
  p.rm_SchnitzNrs=[];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% load SchnitzData
%--------------------------------------------------------------------------
p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
if exist(p.schnitzName)~=2
   error(['Could not read ' p.schnitzName ' , which is required for Schnitz data. Compile first. Exiting ...']);
   return
else
    load(p.schnitzName);
    disp(['Load from ''' p.schnitzName ''' completed...']);
end
if ~existfield(schnitzcells(1),plotField)
  disp(['Field ' plotField ' does not exist. Not plotting!']);
  return;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% If automated colorRange get range from all schnitz data
%--------------------------------------------------------------------------
% create array that contains every data point of plotField for every
% schnitz (including initial schnitzes, sick schnitzes, excluded schnitzes)
if strcmp(p.colorRange,'automatic')==1
    completePlotField=[];
    for i=1:length(schnitzcells)
        idx=find(~isnan(schnitzcells(i).(plotField)));
        if ~isempty(idx)
            completePlotField=[completePlotField, schnitzcells(i).(plotField)(idx)];
        end
    end
    % get mean and std deviation. choose borders by a factor colorWidth*stddev
    % away from mean;
    colorWidth=2;
    stddevPlotField=std(completePlotField);
    meanPlot=mean(completePlotField);
    p.colorRange=[meanPlot-colorWidth*stddevPlotField    meanPlot+colorWidth*stddevPlotField];
    clear completePlotField;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Get existing FluorColor
%--------------------------------------------------------------------------
fluorcolor=p.fluor1;
if p.fluor1=='none' , fluorcolor=p.fluor2;
elseif (p.fluor2=='none' & p.fluor1=='none'), fluorcolor=p.fluor3;
end
expt=genvarname([ 'expt' fluorcolor]);
%--------------------------------------------------------------------------
      

%--------------------------------------------------------------------------
% Guess which frame field corresponds to plotField (either "frames"or
% "Y_frames"
%--------------------------------------------------------------------------

testschnitznr=10;
% plotField exists for each phase image
myframefield=[];
if length(schnitzcells(testschnitznr).frames)==length(schnitzcells(testschnitznr).(plotField))
        myframefield='frames'; % for length, mu, area,...
elseif length(schnitzcells(testschnitznr).frames)+1==length(schnitzcells(testschnitznr).(plotField))
        myframefield='frames'; %necessary for fitted_Y6_mean etc since error in
   %     data association and sometimes 1 entry to long
% plotField exists for Fluor image
else
    fluorframes=[upper(fluorcolor) '_frames'];
    if length(schnitzcells(testschnitznr).(fluorframes))==length(schnitzcells(testschnitznr).(plotField))
        myframefield=fluorframes; % Y6_mean, Y5_mean etc
    elseif length(schnitzcells(testschnitznr).(fluorframes))+1==length(schnitzcells(testschnitznr).(plotField))
        myframefield=fluorframes; % Rates: dY5_sum_dt etc have often one datapoint less
    end
end
if isempty(myframefield)
    disp(['Don''t find appropriate time points (frames) for plotField. I tried fields of schnitz ' num2str(testschnitznr) ...
                '. If bad choice, change ''testschnitznr'' in program. Exiting...']);
    return
else
    disp(['Will use ' myframefield ' as time/frame field.']);
end

%--------------------------------------------------------------------------


 
        
 %--------------------------------------------------------------------------
% loop over every frame number in frameRange
%--------------------------------------------------------------------------
for framenr=p.frameRange
    % load segmentation file
    segfile=[p.segmentationDir, p.movieName, 'seg', str3(framenr), '.mat'];
    if exist(segfile)~=2
        disp(['Could not read ' segfile ' , maybe wrong frame range. Will continue with next frame...']);
    continue
    else
        % prepare for later check if fluorimage exists %Still to check? (NW
        % 2012-03).         
        eval(['clear ' expt]);
        % now really load segfile
        clear Lc LNsub
        load(segfile);
        if ~(exist('Lc') == 1)
            LNsub=Lc; disp(['Frame ', framenr, ' is not approved. Will use LNsub.']);
        end
        % transform to original image size (of phase contrast)
        phaseCropSize=phaseFullSize;
        segImageFull = zeros(phaseCropSize);
        segImageFull((rect(1):rect(3)), (rect(2):rect(4))) = Lc;
        
        % check if frame has fluor image (assumes that all fluor-colours 
        % are taken for the same frame
        eval(['existfluor=exist(''' expt ''');']) % use for later
        
        %--------------------------------------------------------------------------
        % get data from schnitzcell structure
        %--------------------------------------------------------------------------
        % create matrix with 
        % [schnitznr1,cellnr1,cenx_cent1,ceny_cent1,mu1 ; 
        %  schnitznr2,cellnr2,cenx_cent2,ceny_cent2,mu2; 
        %  .... ]
        % cellnr is the cell number in frame framenr. cenx_cent/ceny_cent are x/y
        % positions in a full size image when the centre of mass of the
        % cell colony is centered in the image (see DJK_tracker_djk). mu:
        % growth rate or respectively any other input (e.g. YFP
        % fluorescence). (named "mu", since this will be typical input)
        dataset=[];
        
        %loop over all schnitzes and check if schnitz is in frame "framenr"
        for schnitzrun=1:length(schnitzcells)
            % check if schnitz is excluded
            if ~ismember(schnitzrun,p.rm_SchnitzNrs)
                % if schnitz is not excluded, check if schnitz appears in
                % frame "framnr"
                frameidx=find((framenr+1)==schnitzcells(1,schnitzrun).frames); %+1 corrects for error in frame association (1st entry will e.g. be 44, if cell really appears in 43. shift by 1)
                if ~isempty(frameidx) % schnitzcell exists in frame "framenr"
                    % get data in case of myframefield=frames
                    if strcmp(myframefield,'frames')==1
                        newcellnr=schnitzcells(1,schnitzrun).cellno(frameidx);
                        newx=schnitzcells(1,schnitzrun).cenx_cent(frameidx);
                        newy=schnitzcells(1,schnitzrun).ceny_cent(frameidx);
                        newmu=schnitzcells(1,schnitzrun).(plotField)(frameidx); 
                        dataset=[dataset; schnitzrun, newcellnr, newx, newy, newmu];
                    else % myframefield=Y_frames
                        fluoridx=find((framenr+1)==schnitzcells(1,schnitzrun).(myframefield));
                        % data exists and extra check for last point (there
                        % often Rate data does not exist any more)
                        if (~isempty(fluoridx) & length(schnitzcells(schnitzrun).(plotField))>=fluoridx)
                            newcellnr=schnitzcells(1,schnitzrun).cellno(frameidx); % cellno exists for all frames
                            newx=schnitzcells(1,schnitzrun).cenx_cent(frameidx);  % so do x,y coordinates
                            newy=schnitzcells(1,schnitzrun).ceny_cent(frameidx);
                            newmu=schnitzcells(1,schnitzrun).(plotField)(fluoridx); 
                            dataset=[dataset; schnitzrun, newcellnr, newx, newy, newmu];
                        end
                    end
                end
                     
                  %newcellnr=schnitzcells(1,schnitzrun).cellno(frameidx);
                  %newx=schnitzcells(1,schnitzrun).cenx_cent(frameidx);
                  %newy=schnitzcells(1,schnitzrun).ceny_cent(frameidx);
                  %newmu=schnitzcells(1,schnitzrun).(plotField)(frameidx); 
                  %dataset=[dataset; schnitzrun, newcellnr, newx, newy, newmu];
                
            end
        end
        disp(['Data extraction for frame ' num2str(framenr) ' completed.']);
        %--------------------------------------------------------------------------
        
        %--------------------------------------------------------------------------
        % if data set empty, continue with next frame
        %-------------------------------------------------------------------------- 
        if isempty(dataset)
            disp(['Data for ' plotField ' doesn''t exist for frame ' num2str(framenr) '. Will continue with next frame.']);
            continue
        end
        %--------------------------------------------------------------------------
        
        
        %--------------------------------------------------------------------------
        % calculate distance of each schnitz from centre of colony and
        % correlation (mu,dist)
        %-------------------------------------------------------------------------- 
        centreX_arithm=mean(dataset(:,3));   % arithmetic mean. each schnitz weighed the same, independent of its size
        centreY_arithm=mean(dataset(:,4));
        centreX_mass=round(0.5*size(segImageFull,2));  % centre of mass mean. cenx_cent/ceny_cent are coordinates of schnitzes when centre of mass of colony is right in centre of image
        centreY_mass=round(0.5*size(segImageFull,1));
        
        %dist_arithm=dataset(:,3)-centreX_arithm; % blubb!!! depedence on
        %                                          x-position (y ignored)
        dist_arithm=sqrt((dataset(:,3)-centreX_arithm).^2+(dataset(:,4)-centreY_arithm).^2);
        dist_mass=sqrt((dataset(:,3)-centreX_mass).^2+(dataset(:,4)-centreY_mass).^2);
        mu=dataset(:,5); % redundant. better readability
        
        %******************************************************
        % correlations: adjusted from DJK_plot_scatterColor
        
        % standard deviation (normalized with (length_vector-1))
        sigma_mu=std(mu,1);
        sigma_dist_arithm=std(dist_arithm,1);
        sigma_dist_mass=std(dist_mass,1);
        %disp(['sigma_mu=', num2str(sigma_mu) , '    sigma_dist_arithm=' , num2str(sigma_dist_arithm)]);
        %disp(['cov_arithm=', num2str(cov_arithm)   ,'  corr_mu_dist_arithm=', num2str(corr_arithm)]) ;
        
        % correlation coefficient
        [stat_R_arithm,stat_P_arithm] = corrcoef(mu,dist_arithm);
        stat_R_arithm = stat_R_arithm(2,1);
        stat_P_arithm = stat_P_arithm(2,1);
        [stat_R_mass,stat_P_mass] = corrcoef(mu,dist_mass);
        stat_R_mass = stat_R_mass(2,1);
        stat_P_mass = stat_P_mass(2,1);
        % fitted line BLUBB: also has to be changed to totalleast squares!
        % NW 2012-05
        stat_fitCoef_arithm = polyfit(dist_arithm,mu,1);
        stat_fitCoef_mass = polyfit(dist_mass,mu,1);
        
        %******************************************************
        
        %--------------------------------------------------------------------------
 
        %--------------------------------------------------------------------------
        % 1st plot: segmentation image with colorcoded plotField
        %--------------------------------------------------------------------------
        if show2dimplot==1
            
        % create plotField coded segmentation image, i.e. the matrix entry
        % for each coordinate is the value of plotField (e.g. mu) of the
        % corresponding schnitz
        plotFieldImage=zeros(size(segImageFull));
        %dummy to note where cells are
        blackwhiteCellImage=zeros(size(segImageFull));
        % loop over all cell numbers
        for cellrun=1:max2(segImageFull) % identical to 1:Lc 
            clear whichposition
            whichposition=find(dataset(:,2)==cellrun);
            % necessary for excluded schnitzes
            if ~isempty(whichposition)
                plotFieldImage(segImageFull==cellrun)=dataset(whichposition,5);
                blackwhiteCellImage(segImageFull==cellrun)=1;
            else
                plotFieldImage(segImageFull==cellrun)=0;
            end
        end
        % very crude lifting of low plotField values to above the threshold
        % for the background Color (must be done, because they otherwise
        % appear in background color)
        % necessary because of strange background color construction
        % (TODO in future: might be useful to switch to indexed image)
        numbercolors=size(p.myColorMap,1); DeltaColor=p.colorRange(2)-p.colorRange(1);
        plotFieldImage(plotFieldImage~=0 & plotFieldImage<(p.colorRange(1)+2/numbercolors*DeltaColor))=p.colorRange(1)+2/numbercolors*DeltaColor;
        
        % very crude decreasing of background to a value below zero in case
        % the colormap extends over zero (for production rates). otherwise
        % the background will appear in red or green.
        if p.colorRange(1)<=0
            plotFieldImage(blackwhiteCellImage==0)=p.colorRange(1)-1000;
        end
        
        
        % actual plotting
        fig1=figure();
        set(fig1,'Position', [50 50 size(plotFieldImage,2), size(plotFieldImage,1)]);
        clf
        colormap(p.myColorMap)
        imagesc(plotFieldImage)
        caxis([p.colorRange])
        hold on
        colorbar
        h = colorbar;
        ylabel(h, [plotField],'Interpreter','none','FontSize',13);
        mytitle=[p.movieName '  frame ' num2str(framenr)]; title(mytitle,'FontSize',13);
        % crude check if maybe fluorescence Data which is only interpolated
        % but not of orig fluo image is plotted
        if isempty(strfind(plotField,'mu')) & existfluor==0
            fluowarning=['If this displays Fluorescence Data, then it is only interpolated!!'];
            xlabel(fluowarning,'FontSize',12);
        end
        
        % save image
        figureFileName=[p.movieName '_spatialdependence_' plotField '_frame' str3(framenr)];
        saveSameSize(fig1,'file',[p.NW_saveDir figureFileName '.png'], 'format', 'png');
        disp([' * Saved plot in ' figureFileName '.png']);
        close(fig1);
        end
        %--------------------------------------------------------------------------
        
        
        %--------------------------------------------------------------------------
        % 2nd plot: Scatterplot. correlation between mu and distance to
        % center
        %--------------------------------------------------------------------------
        %
        
       
        % options to format images
        myxlim=[0 max([dist_arithm;dist_mass])*1.02];
        
        % ARITHMETIC MEAN
        if showscatterplot_arithm==1
            
        fig2=figure();
        set(fig2,'Position', [50 50 800 600]);
        clf
        plot(dist_arithm,mu,'.r')
        hold on
        xlim([myxlim])
        xlabel('distance to arithmetic centre of colony','Interpreter','None')
        ylabel([plotField],'Interpreter','None')
        mytitle=[p.movieName '  frame ' num2str(framenr)]; 
        if isempty(strfind(plotField,'mu')) & existfluor==0
            fluowarning=['.  If Fluo, then only interpolated!!'];
            mytitle=[mytitle fluowarning];
        end
        title(mytitle);
        % add fitted line
        stat_fitted_mu = stat_fitCoef_arithm(1)*dist_arithm + stat_fitCoef_arithm(2);
        line( 'Xdata',dist_arithm, ...
                'Ydata',stat_fitted_mu, ...
                'LineStyle','-', ...
                'LineWidth',2, ...
                'Color','k', ...
                'Marker','none');
        % add statistics text
        text(0.02,0.98,['Correlation Coef R / R2 : ' DJK_setDecimalPlaces(stat_R_arithm,2) ' / ' DJK_setDecimalPlaces(stat_R_arithm*stat_R_arithm,2)],'sc' );
        text(0.02,0.94,['Prob R could be 0 (P)   : ' DJK_setDecimalPlaces(stat_P_arithm,4)],                                   'sc');
        text(0.02,0.90,['Regression line         : y = ' DJK_setDecimalPlaces(stat_fitCoef_arithm(1),5) ' * x + ' DJK_setDecimalPlaces(stat_fitCoef_arithm(2),2)], 'sc');
        text(0.40,0.05,['sigma Y: ', num2str(sigma_mu) ,'   noise Y : ' DJK_setDecimalPlaces(sigma_mu/mean(mu),2)], 'sc');
        % save image
        figureFileName=[p.movieName '_spatialcorrelation_arithm_' plotField '_frame' str3(framenr)];
        saveSameSize(fig2,'file',[p.NW_saveDir figureFileName '.png'], 'format', 'png');
        disp([' * Saved plot in ' figureFileName '.png']);
        close(fig2);
        end

        % MASS MEAN
        if showscatterplot_mass==1
        
        fig3=figure();
        set(fig3,'Position', [50 50 800 600]);
        clf
        plot(dist_mass,mu,'.b')
        hold on
        xlim([myxlim])
        xlabel('distance to mass centre of colony','Interpreter','None')
        ylabel([plotField],'Interpreter','None')
        mytitle=[p.movieName '  frame ' num2str(framenr)]; 
        if isempty(strfind(plotField,'mu')) & existfluor==0
            fluowarning=['.  If Fluo, then only interpolated!!'];
            mytitle=[mytitle fluowarning];
        end
        title(mytitle);
        % add fitted line
        stat_fitted_mu = stat_fitCoef_mass(1)*dist_mass + stat_fitCoef_mass(2);
        line( 'Xdata',dist_mass, ...
                'Ydata',stat_fitted_mu, ...
                'LineStyle','-', ...
                'LineWidth',2, ...
                'Color','k', ...
                'Marker','none');
        % add statistics text
        text(0.02,0.98,['Correlation Coef R / R2 : ' DJK_setDecimalPlaces(stat_R_mass,2) ' / ' DJK_setDecimalPlaces(stat_R_mass*stat_R_mass,2)],'sc' );
        text(0.02,0.94,['Prob R could be 0 (P)   : ' DJK_setDecimalPlaces(stat_P_mass,4)],                                   'sc');
        text(0.02,0.90,['Regression line         : y = ' DJK_setDecimalPlaces(stat_fitCoef_mass(1),5) ' * x + ' DJK_setDecimalPlaces(stat_fitCoef_mass(2),2)], 'sc');
        text(0.40,0.05,['sigma Y: ', num2str(sigma_mu) ,'   noise Y : ' DJK_setDecimalPlaces(sigma_mu/mean(mu),2)], 'sc');
        % save image
        figureFileName=[p.movieName '_spatialcorrelation_mass_' plotField '_frame' str3(framenr)];
        saveSameSize(fig3,'file',[p.NW_saveDir figureFileName '.png'], 'format', 'png');
        disp([' * Saved plot in ' figureFileName '.png']);
        close(fig2);        
        end
                
  end
end
        
        
  
%%%%%%%%%%%%% %from Daan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function number_string = DJK_setDecimalPlaces(number, decimal_places);
number_string = sprintf(['%1.' num2str(decimal_places) 'f'], number);

% factor = power(10,decimal_places);
% number = number*factor;
% number = floor(number);
% number = number / factor;
% number_string = num2str(number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
