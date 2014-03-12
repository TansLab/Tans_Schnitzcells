function [schnitzcells] = PN_addToSchnitzes_Phase_at_TimeField_Mu(schnitzcells,phaseField,timeField,newField)
%
% Difference to PN_addToSchnitzes_Phase_at_Timefield.m:  the name of the
% newField is an input argument and can now be chosen manually.
%
% EXAMPLE:
% schnitzcells = PN_addToSchnitzes_Phase_at_TimeField(schnitzcells,['muP15_fitNew_all'],'dY5_time',['muP15_atRate']);
 
% the last argument is the new name of the first argument restricted to the second


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
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
