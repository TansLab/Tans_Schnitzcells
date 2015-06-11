function schnitzcells=NW_correctCellCycle(p,schnitzcells,piecepoly,fieldToCorrect,varargin);
%  Corrects for cell cycle effects (i.e. periodic fluctuations) of data in
%  schnitzcells. Data is corrected with a fitted function piecepoly, which
%  must be obtained for exactly the fieldToCorrect via
%  NW_plot_cellCycleDependence before.
% New schnitzcells is saved in (default):
% /analysisDir/posXcrop-Schnitz.mat
% *********
% correction is performed via subtraction of cell cycle fluctuation and -
% if fluctuation is positive and thus comes from the real data  (Y6_mean
% instead of noise_Y6_mean etc) - also by division
% *******
%
% ************
% Think carefully if you want to use cells that don't have a complete cell
% cycle and can only be corrected by an estimate of next division time.
% Whether it is necessary to take into account these cells seems to depend
% a lot on the actual fitting time for the cross correlations. (i.e. if
% fitTime ends around >=1 cell cycle before experiment-end-time, most of the
% used cells have a complete cell cycle within this range and not many are 
% discarded if we restrict to completeCycle==1
% ************
%
% 
%
% ************
% So far, only the 'real-time-mode' version of NW_plot_CellCycleDependence can
% be corrected (I'll probably not implement the other '0-1-mode', because the one
% here is better.
% ************
%
% OUTPUT
% 'schnitzcells'      schnitzcells with added cycle corrected field
%                     e.g. fieldToCorrect='muP15_fitNew'  ->
%                     'muP15_fitNew_subtr'  'muP15_fitNew_div'
%
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells, where useForPlot should be set
% 'piecepoly'         piecewise polynomial, obtained from
%                     NW_plot_CellCycleDependence
% 'fieldToCorrect'    datafield which shall be corrected. E.g. 'Y6_mean'
%
%
% OPTIONAL ARGUMENTS:
% 'schnitzNameCycCor'  Name under which schnitzcells is saved
% 'restrictToCompleteCycle'  0: for schnitzes without complete cell cycle
%                    an estimate of its future division time and thus a
%                    correction for cell cycle dependence is made with help
%                    of av_mu_fitNew and average increase factor of length
%                    (typically =2). If estimated time is too far off mean,
%                    then a threshold condition restricts it to
%                    X*mean_interdivision_time
%                    1: only cells with complete cell cycle are used
%                    default: 1
% NW_saveDir         where corrected schnitzcells structure should be
%                    stored.
%                    default: posXcrop/data/
% 'mytimefield'      timefield corresponding to fieldY (must have sdame length).
%                    E.g. 'time' or 'time_atdR'. Default: 'time'
%
%


%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'NW_correctCellCycle'; 

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
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(p,'schnitzNameCycCor')
 % p.schnitzNameCycCor = [p.tracksDir,p.movieName,'-Schnitz-CycCor.mat'];
 p.schnitzNameCycCor = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end
if ~existfield(p,'restrictToCompleteCycle')
  p.restrictToCompleteCycle = 1;
end

% Just in case it has not been set yet
if ~existfield(p,'NW_saveDir')
  p.NW_saveDir = [p.tracksDir];
end
% Make sure that DJK_saveDir directory exists
if exist(p.NW_saveDir)~=7
  [status,msg,id] = mymkdir([p.NW_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.NW_saveDir ' : ' msg]);
    return;
  end
end

if ~existfield(p,'mytimefield')    
  p.mytimefield = 'time';
end
mytimefield=p.mytimefield; % easier to write down
% check if total length of datafield and timefield are identical
alldata=[schnitzcells.(fieldToCorrect)];
if strcmp('time',mytimefield)~=1
    idx=find(~isnan(alldata)); % muP15_fitNew can contain 'NaN' data but no correpsonding time point data
                               % for muP15_fitNew_all, the NaN values have
                               % to stay... argh
    alldata=alldata(idx);
end
alltime=[schnitzcells.(mytimefield)];
if length(alldata)~=length(alltime)
    disp(['Datafield and timefield don''t correspond (have different size).'])
    return
end
clear alldata alltime
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Can''t do anything!');
  return;
end

if ~existfield(schnitzcells(1),fieldToCorrect)
  disp(['Field ' fieldToCorrect ' does not exist. Not plotting!']);
  return;
end

% check if field "frames" for time calculation exists
if (~existfield(schnitzcells(1),'frames'))
  disp(['Field ''frames'' does not exist. Cannot calculate Time. Exiting...']);
  return;
end

% if useForPlot has not been set
useAllcells = ~existfield(schnitzcells(1),'useForPlot'); %if useForPlot has not been set at all

% if cells with non-complete cycle are corrected, we need the fields "time" ,
% 'av_birthLength_fitNew','birthTime', 'av_divLength_fitNew','interDivTime'
% "av_mu_fitNew" field (the may be 'NaN')
if p.restrictToCompleteCycle==0
    if (~existfield(schnitzcells(1),'av_birthLength_fitNew'))
      disp(['Field ''av_birthLength_fitNew'' does not exist. Exiting...']);
      return;
    end
    if (~existfield(schnitzcells(1),'birthTime'))
      disp(['Field ''birthTime'' does not exist. Exiting...']);
      return;
    end
    if (~existfield(schnitzcells(1),'av_divLength_fitNew'))
      disp(['Field ''av_divLength_fitNew'' does not exist. Exiting...']);
      return;
    end
    if (~existfield(schnitzcells(1),'interDivTime'))
      disp(['Field ''interDivTime'' does not exist. Exiting...']);
      return;
    end
    if (~existfield(schnitzcells(1),'av_mu_fitNew'))
      disp(['Field ''av_mu_fitNew'' does not exist. Exiting...']);
      return;
    end
    if (~existfield(schnitzcells(1),'time'))
      disp(['Field ''time'' does not exist. Exiting...']);
      return;
    end
end
%--------------------------------------------------------------------------

% %--------------------------------------------------------------------------
% % get time field which corresponds to fieldToCorrect
% %--------------------------------------------------------------------------
% % if fieldToCorrect exists for all frames, take "frames" as timefield
% testschnitznr=100;
%if length(schnitzcells(testschnitznr).frames)==length(schnitzcells(testschnitznr).(fieldToCorrect))
%     mytimefield='frames';
   
        
% % if fieldToCorrect exists for less frames, try with fluo-frames as timefield
%else
%     %get 1st real fluor timefield

%if strcmp(p.fluor1,'none')==0
 %         mytimefield=[upper(p.fluor1) '_frames'];
 %    elseif strcmp(p.fluor2,'none')==0
 %         mytimefield=[upper(p.fluor2) '_frames'];
 %    elseif strcmp(p.fluor3,'none')==0
 %         mytimefield=[upper(p.fluor2) '_frames'];
 %    else
 %         disp(['Don''t find appropriate time points (frames) for fieldToCorrect. No FluorColor? Exiting...']);
 %         return
%     end
%     if length(schnitzcells(testschnitznr).(mytimefield))==length(schnitzcells(testschnitznr).(fieldToCorrect))
%          disp(['Will use ' mytimefield ' as time/frame field.'])
%     else
%          disp(['Don''t find appropriate time points (frames) for fieldToCorrect. Exiting...']);
%          disp(['I checked fields of schnitz ' num2str(testschnitznr) '. Maybe bad choice. Change in code testschnitznr=...'])
%          return
%     end
%end
% %--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% if also using cells without complete cell cycle, get average length
% increase factor per complete cell cycle
%--------------------------------------------------------------------------
if p.restrictToCompleteCycle==0
    divlength=[];   % array with division lengths (of schntizes with complete cell cycle)
    birthlength=[]; % array with birth lengths
    mu=[];          % array with average overall growth rate, fitted to length increase of complete cycle (av_mu_fitNew)
    lifetime=[];    % array with life time of each schnitz
    numberframes=[];% array with number of frames a cell appears in
    
    for schnitz=1:length(schnitzcells)
        if (useAllcells | schnitzcells(schnitz).useForPlot==1) & schnitzcells(schnitz).completeCycle==1
            % should always be the case ...
            if ~isnan(schnitzcells(schnitz).av_birthLength_fitNew) & ~isnan(schnitzcells(schnitz).av_divLength_fitNew) %should always be the case
                birthlength=[birthlength; schnitzcells(schnitz).av_birthLength_fitNew];
                divlength=[divlength; schnitzcells(schnitz).av_divLength_fitNew];
                mu=[mu; schnitzcells(schnitz).av_mu_fitNew];
                lifetime=[lifetime; schnitzcells(schnitz).interDivTime];
                numberframes=[numberframes;length(schnitzcells(schnitz).frames)];
            %else
            %    disp(schnitz)
            end
        end
    end

    lengthfrac=divlength./birthlength; %factor for length increase of each schnitz
    % plot figure
    showhistfig=0;
    if showhistfig==1
        figure
        clf
        hist(numberframes) %hist(lengthfrac)
        axis on
    end
        
    meanlengthfrac=mean(lengthfrac); %average length factor increase for all schnitzes
                               % should be very close to 2
    meanbirth=mean(birthlength);
    meandiv=mean(divlength);
    meanmu=mean(mu);
    meantime=mean(lifetime);
    meannumberframes=mean(numberframes);
    disp(['mean length increase: factor ' num2str(meanlengthfrac)])
    disp(['mean birth length: ' num2str(meanbirth) ', mean division length: ' num2str(meandiv)])
    disp(['mean growth rate (av_mu_fitNew) ' num2str(meanmu) ', mean life time [min] ' num2str(meantime)]);
    disp(['mean number of frames per schnitz: ' num2str(meannumberframes)])
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% create new field names and prepare correction data
%--------------------------------------------------------------------------
% new field names
field_subtr=[fieldToCorrect, '_subtr']; %correction subtracted
field_div=[fieldToCorrect, '_div']; %correction divided 
%get average value of spline interpolation of fieldToCorrect 
% -> normalize
xx=0:0.01:1;
yy=ppval(piecepoly,xx);


% subtract cell cycle noise (force mean(noise) to be=0)
SplineMeanToSubstract=mean(yy);

% If polynomial is always above zero, it results from real Data (not
% noise_XY) and a _div field can be added.
% divide by cell cycle noise (mean(noise) forced to be =1)
if ~isempty(find(yy<0))
    disp(['Cannot add _div division corrected field for ' fieldToCorrect]);
    addDivision=0;
else 
    disp(['Will also add _div division corrected field for ' fieldToCorrect]);
    addDivision=1;
    SplineMeanToDivide=SplineMeanToSubstract; % redundant...
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% loop over all schnitzes and create corrected field
%--------------------------------------------------------------------------
% ALWAYS CALCULATES TIME IN 'real-time-mode' of
% NW_plot_CellCycleDependence!
% '0-1-mode' not implemented

for schnitz=1:length(schnitzcells) 
    if useAllcells | schnitzcells(schnitz).useForPlot % test: use schnitz
        
      % ************************
        if ~isempty(schnitzcells(schnitz).(fieldToCorrect)) % sometimes weird (long) cells don't have fluo frames
            
            % ------------------------------------------------
            % CORRECT SCHNITZES WITH COMPLETE CELL CYCLE
            % ------------------------------------------------
            if schnitzcells(schnitz).completeCycle==1     
                clear dataY totalTime YdataTime dataY_subtr dataY_div

           %     % number of frames a schnitz exists (frames used as proxy for time)    
           %     totalTime=length(schnitzcells(schnitz).frames); 
           %     % absolute time points (in frames) of fieldToCorrect with initial 
           %     % birth frame set to 0
           %     YdataTime=schnitzcells(schnitz).(mytimefield)-schnitzcells(schnitz).frames(1); 
           %     % relative time (birth=0. last frame=(#frames-1)/#frames)
           %     YdataTime=YdataTime/totalTime; 

                 totalTime=schnitzcells(schnitz).divTime-schnitzcells(schnitz).birthTime;
                YdataTime=schnitzcells(schnitz).(mytimefield)-schnitzcells(schnitz).birthTime; %time points in [min] with respect to birth time of cell
                YdataTime=YdataTime/totalTime; % time points normalized to interval [0 1]  -> phase
        
         
                
                %get fieldToCorrect data of this schnitz
                dataY=schnitzcells(schnitz).(fieldToCorrect);


                %create arrays for corrected data
                dataY_subtr=[]; %zeros(length(dataY),1);
                if addDivision==1
                    dataY_div=[]; %zeros(length(dataY),1);
                end

          
   
                if length(YdataTime)>0  %should always be the case for complete cycle
       %             % STILL TO TEST! BLUBB!
       %             % If dataY is one shorter than YdataTime, remove
       %             % last time point of YdataTime, since last (daughter) datapoint of dataY is missing 
       %             % (happens when: 
       %             % - not complete cell cycle
       %             % - complete celly cycle, but daughters don't have
       %             % any fluor image
       %             if length(YdataTime)==length(dataY)+1
       %                   YdataTime(end)=[];
       %             end       
       %             % ****

                    % correction by subtraction
                    dataY_subtr=dataY-ppval(piecepoly,YdataTime)+SplineMeanToSubstract;
                    % % correction by division
                    %dataY_div=dataY./(ppval(piecepoly,YdataTime)+forceShift).*SplineMeanToDivide;  %DOES NOT WORK YET
                    if addDivision==1
                        dataY_div=dataY./ppval(piecepoly,YdataTime).*SplineMeanToDivide; 
                    end
                end

                % add to schnitzcells
                schnitzcells(schnitz).(field_subtr)=dataY_subtr;
                if addDivision==1
                        schnitzcells(schnitz).(field_div)=dataY_div;
                end
                 
            % ------------------------------------------
            % TREAT SCHNITZES WITH INCOMPLETE CELL CYCLE   
            % ------------------------------------------
            % here, the 'div' case is ignored
            elseif schnitzcells(schnitz).completeCycle==0    
                % only consider complete cycle -> associate empty array
                % (useful?)
                % -----------------------
                % don't add data
                % -----------------------
                if p.restrictToCompleteCycle==1
                    schnitzcells(schnitz).(field_subtr)=[]; % if problems: change to array of NaN (NW 2012-03)

                % ----------------------
                % correct incomplete schnitzcell cycles by estimate
                % -----------------------
                elseif p.restrictToCompleteCycle==0
                    clear dataY totalTime YdataTime dataY_subtr dataY_div

                    % ------------------------------
                    % no average growth rate exists for cell
                    % ------------------------------
                    if isnan(schnitzcells(schnitz).av_mu_fitNew) %very late cells (or first cell)  ->no growth rate
              %          %partiallifetime=1/meannumberframes; 
              %          %set its growth rate to average
              %          %for this case, unit is NOT in minutes!
              %          partiallifetime=length(schnitzcells(schnitz).frames)+1;
              %          predictedlifetime=meannumberframes;
                         if ~isnan(schnitzcells(schnitz).birthTime)
                             partiallifetime=schnitzcells(schnitz).time(end)-schnitzcells(schnitz).birthTime;
                         else
                             partiallifetime=schnitzcells(schnitz).divTime-schnitzcells(schnitz).time(1);
                         end
                        predictedlifetime=meantime;
                        totalTime=meantime;
                        disp(['schnitz ' num2str(schnitz) ': measured time=' num2str(partiallifetime) ...
                                 ' estimated life time=' num2str(meantime)]);
                    
                    % ---------------------------
                    % average growth rate exists
                    % ---------------------------
                    else  
                        % time [min] since schnitz was born: partiallifetime

                        % SET BIRTHTIME OF FIRST SCHNITZ TO=0
                        % ROUGH ESTIMATE BUT SCHNITZ WILL PROBABLY NEVER BE
                        % USED AND OTHERWISE DJK_getBranches collapses
                        if isnan(schnitzcells(schnitz).birthTime)
                            partiallifetime=schnitzcells(schnitz).time(end); %...
                                %+ 0.5*(schnitzcells(schnitz).time(2)-...
                            %schnitzcells(schnitz).time(1));
                        else
                            partiallifetime=schnitzcells(schnitz).time(end)- schnitzcells(schnitz).birthTime ;
                        %+ 0.5*(schnitzcells(schnitz).time(2)-schnitzcells(schnitz).time(1));
                         
                        end

                        % time, after which cell theoretically would have length
                        % increase by factor "meanlengthfrac"
                        % meanlengthfrac = 2^(lifetime/60*mu)
                        % **** Here, the estimate comes in, that every cell is
                        % supposed to grow by a factor meanlengthfrac before
                        % division ******
                        predictedlifetime=60*log(meanlengthfrac)/(log(2)*schnitzcells(schnitz).av_mu_fitNew);
                    
                        % in case cell is already longer, assumes that it's gonna
                        % divide in the next frame
                        if predictedlifetime<=partiallifetime
                            % number of frames in whic schnitz is supposed to exist
               %             totalTime=length(schnitzcells(schnitz).frames); 
                            totalTime=partiallifetime; 

               %             %tell user
                            disp(['schnitz ' num2str(schnitz) ': measured time=' num2str(partiallifetime) ...
                                 ' estimated life time= id']);
                        else
                            totalTime=predictedlifetime;

                            % if the #existing time points is too low, the estimate
                            % can be very crappy.
                            % if less then 50% of average life time recorded
                            % and estimated lifetime higher/lower than
                            % 1.3*max or 1*min (cells get slower in end of experiment!) of
                            % complete cell cycle, then reset value to these border values
                            if (partiallifetime<0.5*meantime) & (totalTime>1.3*max(lifetime))
                                    totalTime=1.3*max(lifetime);
                            elseif (partiallifetime<0.5*meantime) & (totalTime<min(lifetime))
                                    totalTime=min(lifetime);
                            end                    

                            %tell user
                            disp(['schnitz ' num2str(schnitz) ': measured time=' num2str(partiallifetime) ...
                                 ' estimated life time=' num2str(totalTime)]);
                        end % end compare predicted lifetime & partiallifetime
                    end % disinguish if av_mu_fitNew exists

                        
             %           % position in cell cycle (0 to 1)
             %           phi=partiallifetime/predictedlifetime;
             %           % # frames corresponding to phi
             %           partialframes=length(schnitzcells(schnitz).frames)-1;
             %           % estimate number of frames cell will live in
             %           allframesestimate=round(partialframes/phi);

             %           % get proxy for total Time
             %           totalTime=1+allframesestimate;

              %          % if the #existing frames is too low, the estimate
              %          % can be very crappy.
              %          % if less then 50% of average life time recorded
              %          % and estimatedframenumber higher/lower than
              %          % 1.3*max or 0.5*min (cells get slower in end!) of
              %          % complete cell cycle, then stick to these border
              %          % values
              %          if ((partialframes+1)<0.5*meannumberframes) & (totalTime>1.3*max(numberframes))
              %                  totalTime=1.3*max(numberframes);
              %          elseif ((partialframes+1)<0.5*meannumberframes) & (totalTime<min(numberframes))
              %                  totalTime=min(numberframes);
              %          end                    

              %          %tell user
              %          disp(['schnitz ' num2str(schnitz) ': actual frames=' num2str(partialframes+1) ...
              %               ' estimated frames=' num2str(totalTime)])
              %      end

              % -------------------------------------------------------------------------
                        % *****************************************
                        % now proceed like in completeCycle==1 case
                        % *****************************************
              % --------------------------------------------------------------------
             %           % absolute time points (in frames) 
             %           YdataTime=schnitzcells(schnitz).(mytimefield)-schnitzcells(schnitz).frames(1); 
             %           % relative time (birth=0. last frame=(#frames-1)/#frames)
             %           YdataTime=YdataTime/totalTime; 
                        
                       % absolute time points (in frames)
                       if ~isnan(schnitzcells(schnitz).birthTime)
                           YdataTime=schnitzcells(schnitz).(mytimefield)-schnitzcells(schnitz).birthTime;
                       else
                           YdataTime=schnitzcells(schnitz).(mytimefield); % ROUGH ESTIMATE: SCHNITZ 1 IS BORN AT TIME=0
                       end
                        % relative time (birth=0. last frame=(#frames-1)/#frames)
                        YdataTime=YdataTime/totalTime; 

                        %get fieldToCorrect data of this schnitz
                        dataY=schnitzcells(schnitz).(fieldToCorrect);

                        %create arrays for corrected data
                        dataY_subtr=[]; %zeros(length(dataY),1);
                        if addDivision==1
                            dataY_div=[]; %zeros(length(dataY),1);
                        end

                %        if length(YdataTime)>0 %should always be the case for complete cycle
                %            % STILL TO TEST! BLUBB!
                %            % If dataY is one shorter than YdataTime, remove
                %            % last time point of YdataTime, since last (daughter) datapoint of dataY is missing 
                %            % (happens when: 
                %            % - not complete cell cycle
                %            % - complete celly cycle, but daughters don't have
                %            % any fluor image
                %            if length(YdataTime)==length(dataY)+1
                %                YdataTime(end)=[];
                %            end
                %            % ***

                            if ~isempty(dataY) & ~isnan(dataY) % isnan check: muP15_fitNew can exist as 'nan' vector even though no fluor time point exists
                                % correction by subtraction
                                dataY_subtr=dataY-ppval(piecepoly,YdataTime)+SplineMeanToSubstract;
                                if addDivision==1 %correct for division
                                    dataY_div=dataY./ppval(piecepoly,YdataTime).*SplineMeanToDivide; 
                                end
                            end
                end % end of: restrictToCompleteCycle yes/no

                        % add to schnitzcells
                        schnitzcells(schnitz).(field_subtr)=dataY_subtr;
                        if addDivision==1
                            schnitzcells(schnitz).(field_div)=dataY_div;
                        end
            end
     
        else % only for weird long cells
           schnitzcells(schnitz).(field_subtr)=[];
           if addDivision==1
                  schnitzcells(schnitz).(field_div)=[];
           end 
        end % end weird long cells *****

    end
end


%--------------------------------------------------------------------------
% Save new schnitzcells with corrected data added 
%--------------------------------------------------------------------------
save(p.schnitzNameCycCor, 'schnitzcells');
disp('');
disp(['Save in ''' p.schnitzNameCycCor ''' completed...']);
%--------------------------------------------------------------------------




