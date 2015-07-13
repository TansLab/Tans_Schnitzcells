%
%schnitz= toschnitz(cellid,schnitz,fieldname, data)
%
%adds a field '<fieldname>' into the schnitz cell-number-labeled structure array
%from the data indexed in the format: {frame_number}(label-number)
%cellid- label-cell number correspondence array
%T.B. 04/05
%
% This function has two modes, if P is only one parameter, it will be
% assumed that each of its values corresponds to a schnitz.
% If P is a cell, and contains multiple parameters, it is assumed the input
% parameter cellid can be used to link a frame number to cell ids.
% - Comment MW

function   schnitz= toschnitz(cellid,schnitz,f,P)

% If P contains multiple fields
if(iscell(P))

    %% Initializing
    % create schnitz if doesn't exist
    if(~isstruct(schnitz))
        schnitz=[]
    end

    % Name dummy if no name given
    if(isempty(f))
        f='dummy';
    end

    % loop over schnitzes and add arrays w. fieldnames
    for i=1:size(schnitz,2)
        schnitz(i).(f)=[];
    end

    %% now for each parameter in P 
    for i=1:size(P,2)

        for k=1:numel(P{i})

            if( i > size(cellid,2) | k > max(size(cellid{i})) )
                warning(sprintf('input array index {%d}(%d) exceeds cellid array size, skipped',i,k))
            else                
                id=cellid{i}(k);                                

                if(id>size(schnitz,2))

                    schnitz(id).(f)=P{i}(k);

                else
                    schnitz(id).(f)=[schnitz(id).(f) P{i}(k)];
                end
            end
        end
    end

    if(strcmp(f,'dummy'))

        fn=fieldnames(schnitz(1).dummy);
        for i=1:size(schnitz,2)

            if isempty(schnitz(i).dummy)
               continue
            end

            for j=1:size(fn,1)
                schnitz(i).(fn{j})=[schnitz(i).dummy.(fn{j})];
            end
        end

        schnitz=rmfield(schnitz,'dummy');
    end

else % if P contains single field
    
    %% each schnitz gets the field P, without using lookup table
    for i=1:max(size(P))
        schnitz(i).(f)=P(i);
    end

end