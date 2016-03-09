function [schnitzcells] = PN_addToSchnitzes_Phase_at_TimeField(schnitzcells,phaseField,timeField)


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
newField = [phaseField '_at_' timeField];
for i = 1:length(schnitzcells)
    schnitzcells(i).(newField) = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);

    if ~isempty(s.(timeField))
        for f = 1:length(s.(timeField))
            age = s.time == s.(timeField)(f);
            s.(newField)(end+1) =  s.(phaseField)(age);
        end
    end

    schnitzcells(i) = s;
end
