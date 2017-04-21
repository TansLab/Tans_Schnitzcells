function [myTraces, myDivisions, traceIdxs] = MW_gettracefromschnitzcells(schnitzcells,STARTSCHNITZES,TIMEFIELD,YFIELD,criterium)
% criterium should be longest, shortest, random, or vector giving the
% daughters

traceIdxs=[];
for startSchnitzIdx = STARTSCHNITZES
        
    % initialize trace
    myTraces{startSchnitzIdx} = [];
    myDivisions{startSchnitzIdx} = [];
    % initialize loop
    timeIdx = 1;
    schnitzIdx = startSchnitzIdx;

    counter=1;
    while schnitzIdx ~= 0

        % if there was a division
        if timeIdx > numel(schnitzcells(schnitzIdx).(YFIELD))

            % get two daughters
            daughter1Idx=schnitzcells(schnitzIdx).D;
            daughter2Idx=schnitzcells(schnitzIdx).E;

            % if end of lineage quit making trace
            if daughter1Idx==0 && daughter2Idx==0
                break;
            end                                    

            % determine the shortest/longest daughter
            if schnitzcells(daughter1Idx).(YFIELD)(1)>schnitzcells(daughter2Idx).(YFIELD)(1)
                longestSchnitzIdx = daughter1Idx;
                shortestSchnitzIdx = daughter2Idx;
            else
                longestSchnitzIdx = daughter2Idx;
                shortestSchnitzIdx = daughter1Idx;
            end
            
            % select the daughter based on criterium
            if isnumeric(criterium) % specific trace pre-determined
                schnitzIdx = criterium(counter);
                counter=counter+1;
            elseif strcmp(criterium,'longest')
                schnitzIdx = longestSchnitzIdx; 
            elseif strcmp(criterium,'shortest'), schnitzIdx = shortestSchnitzIdx; 
            elseif strcmp(criterium,'random')
                toChooseFrom=[daughter1Idx,daughter2Idx];
                schnitzIdx  = toChooseFrom(double(rand()<0.5)+1);
                traceIdxs(end+1) = schnitzIdx;
            end
            
            % log the division event in a "trace" of its own
            myDivisions{startSchnitzIdx} = [myDivisions{startSchnitzIdx}; schnitzcells(schnitzIdx).(TIMEFIELD)(1), schnitzcells(schnitzIdx).(YFIELD)(1)];        

            % reset time counter for schnitz
            timeIdx = 1;

        end

        % create trace            
        myTraces{startSchnitzIdx} = [myTraces{startSchnitzIdx}; schnitzcells(schnitzIdx).(TIMEFIELD)(timeIdx), schnitzcells(schnitzIdx).(YFIELD)(timeIdx)];        

        % increase time counter for schnitz
        timeIdx = timeIdx +1;                            

    end                

end