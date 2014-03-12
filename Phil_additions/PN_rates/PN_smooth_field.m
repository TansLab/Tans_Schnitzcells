% DJK_addToSchnitzes_atYandDY
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'field'             name of field
% 'varType'           'intensive' or 'extensive'. Determines how parent and
%                     daughter information is taken into account
%                     (divided/multiplied by 2 or not)
%                     default: 'intensive'
%
function [schnitzcells] = PN_smooth_field(schnitzcells,field,varType)

% check input (very basic)
if strcmp(upper(varType),'INTENSIVE')==1
    typeInt=true;
elseif strcmp(upper(varType),'EXTENSIVE')==1
    typeInt=false;
else
    disp('Don''t recognize variable type (intensive, extensive). Will assume intensive.')
    typeInt=true;
end


%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
smooth_field                     = [field '_s'];
for i = 1:length(schnitzcells)
    schnitzcells(i).(smooth_field)   = [];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    field_temp = s.(field);
    
    if ~isempty(field_temp)
        beginning_time = 1;
        ending_time = length(field_temp);
        clear sD sE sP;
        
        %add parent
        if s.P ~= 0
            sP = schnitzcells(s.P);
            if ~isempty(sP.(field))
                if typeInt %intensive
                    field_temp = [sP.(field)(end) field_temp];
                else %extensive
                    field_temp = [sP.(field)(end)/2 field_temp]; %assumes equal distribution btw daughters
                end
                beginning_time = 2;
            end
        end
        
        %add children
        if s.D ~= 0 && s.E ~= 0
            sD = schnitzcells(s.D);
            sE = schnitzcells(s.E);
            if ~isempty(sD.(field)) && ~isempty(sE.(field))
                if typeInt
                    field_temp(end+1) = (sD.(field)(1) + sE.(field)(1))/2;
                else
                    field_temp(end+1) = sD.(field)(1) + sE.(field)(1);
                end
            end
        end
        
        field_temp = smooth(field_temp,'lowess');
        ending_time = ending_time + beginning_time-1;
        field_temp = field_temp(beginning_time:ending_time);
    end
    
    s.(smooth_field) = field_temp';
    schnitzcells(i) = s;
end












