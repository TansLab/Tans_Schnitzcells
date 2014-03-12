%
%schnitz= toschnitz(cellid,schnitz,fieldname, data)
%
%adds a field '<fieldname>' into the schnitz cell-number-labeled structure array
%from the data indexed in the format: {frame_number}(label-number)
%cellid- label-cell number correspondence array
%T.B. 04/05

function   schnitz= toschnitz(cellid,schnitz,f,P)

if(iscell(P))


    if(~isstruct(schnitz))
        schnitz=[]
    end

    if(isempty(f))
        f='dummy';
    end

    for i=1:size(schnitz,2)
        schnitz(i).(f)=[];
    end

    for i=1:size(P,2)
        for k=1:max(size(P{i}))


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

else

    for i=1:max(size(P))
        schnitz(i).(f)=P(i);
    end

end