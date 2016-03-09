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
% 'lengthField'       sets which field of schnitzcells to use for length
%                     determination, usually 'length_fitNew'
%
function [schnitzcells] = PN_smooth_field(schnitzcells,field,varType,lengthField)

% check input (very basic)
if strcmp(upper(varType),'INTENSIVE')==1
    typeIntensive=true;
elseif strcmp(upper(varType),'EXTENSIVE')==1
    typeIntensive=false;
else
    disp('Don''t recognize variable type (intensive, extensive). Will assume intensive.')
    typeIntensive=true;
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
                if typeIntensive %intensive
                    field_temp = [sP.(field)(end) field_temp];
                else %extensive
                    
                    % get ratio of length daughter:parent. -MW 2015-07
                    thisGuysLengths = s.(lengthField);
                    parentLengths = sP.(lengthField);
                    lengthRatio = thisGuysLengths(1)/parentLengths(end);
                    
                    field_temp = [sP.(field)(end)*lengthRatio field_temp]; %assumes equal distribution btw daughters
                end
                beginning_time = 2;
            end
        end
        
        %add children
        if s.D ~= 0 && s.E ~= 0
            sD = schnitzcells(s.D);
            sE = schnitzcells(s.E);
            if ~isempty(sD.(field)) && ~isempty(sE.(field))
                if typeIntensive
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












