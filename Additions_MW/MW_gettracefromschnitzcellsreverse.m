function [myTraces] = MW_gettracefromschnitzcellsreverse(schnitzcells,lastSchnitzes,timeField,yFields)
% function [myTraces] = MW_gettracefromschnitzcellsreverse(schnitzcells,lastSchnitzes,timeField,yFields)
%
% Input:
%   schnitzcells    standard schnitzcells struct
%   lastSchnitzes   the schnitzes with which your traces of interest end
%                   e.g. [100,110,150]
%                   Leave empty, i.e. [], to automatically take all ends
%   timeField       the field that should be used for time data
%                   e.g. 'time'
%   yFields         the fields of interest for which you want values
%                   e.g. {'muP9_fitNew_all','Y5_mean_all'}
%
%
% Output:
% 
%     myTraces(i).lineageSchnitzNrs     schnitznr's of mother, daugther, granddaughter, grand-granddaughter etc.
%     myTraces(i).time                  your timeField
%     myTraces(i).divisionIndices       indices of the trace at which divisions occur
%     myTraces(i).('YOURFIELD')         values of the fields you supplied in yFields


%% Initialize

% Take all last schnitzes if none are given
if isempty(lastSchnitzes)
   lastSchnitzes = find([schnitzcells.D]==0);
end

%% Create schnitz lineage tree
myTraces=struct;
traceCounter = 0;
for lastSchnitzIdx = lastSchnitzes
 
    %% Initialize
    traceCounter=traceCounter+1;
    schnitzIdx = lastSchnitzIdx;
    myTraces(traceCounter).lineageSchnitzNrs = [];
    
    
    counter=1;
    while schnitzIdx ~= 0
        
        % store current schnitz
        myTraces(traceCounter).lineageSchnitzNrs(end+1) = schnitzIdx;

        % get parent
        schnitzIdx=schnitzcells(schnitzIdx).P;

        % if end of lineage quit making trace
        if schnitzIdx==0 
            break;
        end                                    

    end
       
    %% Create the chronological tree
    chronologicalTree = myTraces(traceCounter).lineageSchnitzNrs;
    chronologicalTree = chronologicalTree(end:-1:1);
    myTraces(traceCounter).lineageSchnitzNrs=chronologicalTree;
    
    %% Compile traces
    myTraces(traceCounter).(timeField)=[];    
    myTraces(traceCounter).divisionIndices=[];
    for fIdx = 1:numel(yFields)
        yFieldCurrent=yFields{fIdx};
        myTraces(traceCounter).(yFieldCurrent)=[];
    end
    for schnitzIdx = myTraces(traceCounter).lineageSchnitzNrs
        
        % TIMEFIELD
        myTraces(traceCounter).(timeField) = [myTraces(traceCounter).(timeField), schnitzcells(schnitzIdx).(timeField)];
        
        % YFIELDS
        if ~isempty(schnitzcells(schnitzIdx).(yFieldCurrent))
            % for each y-field of interest create a trace
            for fIdx = 1:numel(yFields)
                yFieldCurrent=yFields{fIdx};
                myTraces(traceCounter).(yFieldCurrent)    = [myTraces(traceCounter).(yFieldCurrent), schnitzcells(schnitzIdx).(yFieldCurrent)];        
            end            
        else
            % if there is no data, create nan values for missing data
            % (E.g. a 1-timepoint schnitzcell has no growth-rate
            % calculated, but does give timepoints
            for fIdx = 1:numel(yFields)
                yFieldCurrent=yFields{fIdx};
                myTraces(traceCounter).(yFieldCurrent)    = [myTraces(traceCounter).(yFieldCurrent), nan(1,numel(schnitzcells(schnitzIdx).(yFieldCurrent)))];        
            end            
        end
        
        % DIVISIONS
        myTraces(traceCounter).divisionIndices = [myTraces(traceCounter).divisionIndices, zeros(1,numel(schnitzcells(schnitzIdx).(timeField))-1) 1];
        
    end
        
    
end
    

%%

%{
%%
myTraces=struct;
traceCounter = 1;
for lastSchnitzIdx = lastSchnitzes
        
    % initialize trace
    myTraces(traceCounter).(YFIELD) = [];
    myTraces(traceCounter).(TIMEFIELD) = [];
    %myTraces(traceCounter).divisionTimepoints = [];
    myTraces(traceCounter).endSchnitz = lastSchnitzIdx;
    myTraces(traceCounter).divisionIndices = [];
    myTraces(traceCounter).lineageSchnitzNrsInverse = [];
 
    % initialize loop    
    schnitzIdx = lastSchnitzIdx;
    %timeIdx = numel(schnitzcells(schnitzIdx).(TIMEFIELD)(end));

    counter=1;
    while schnitzIdx ~= 0

        % if there was a division
        if timeIdx < 1

            % store current schnitz
            myTraces(traceCounter).lineageSchnitzNrsInverse = schnitzIdx;
            
            % get parent
            schnitzIdx=schnitzcells(schnitzIdx).P;

            % if end of lineage quit making trace
            if schnitzIdx==0 
                break;
            end                                    

            % log the division event in a "trace" of its own            
            %myTraces(traceCounter).divisionTimepoints = [myTraces(traceCounter).divisionTimepoints; schnitzcells(schnitzIdx).(TIMEFIELD)(end), schnitzcells(schnitzIdx).(YFIELD)(1)];
            
            % log the indices of division events
            myTraces(traceCounter).divisionIndices = numel(myTraces(traceCounter).(TIMEFIELD))+1;
            
            % reset time counter for schnitz
            %timeIdx = numel(schnitzcells(schnitzIdx).(TIMEFIELD)(end));

        end

        % create trace            
        myTraces(traceCounter).(TIMEFIELD) = [myTraces(traceCounter).(TIMEFIELD), schnitzcells(schnitzIdx).(TIMEFIELD)];
        if ~isempty(schnitzcells(schnitzIdx).(YFIELD))
            myTraces(traceCounter).(YFIELD)    = [myTraces(traceCounter).(YFIELD), schnitzcells(schnitzIdx).(YFIELD)];        
        else
            myTraces(traceCounter).(YFIELD)    = [myTraces(traceCounter).(YFIELD), nan(1,numel(schnitzcells(schnitzIdx).(YFIELD)))];        
        end

        % increase time counter for schnitz
        %timeIdx = timeIdx -1;                            

    end                

    traceCounter=traceCounter+1;
    
end

%{
IDX=1;

figure; clf; hold on;
plot(myTraces(IDX).time,myTraces(IDX).muP9_fitNew_all,'-o','LineWidth',2);


%}
%}




