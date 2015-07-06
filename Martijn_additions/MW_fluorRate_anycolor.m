function [schnitzcells] = MW_fluorRate_anycolor(schnitzcells, fluorcolor, lengthField)
% Adds rate for any color
% DJK_addToSchnitzes_fluorRate_phase.
% 
% input
%   - color    
%       letter that represents fluors color
%   - lengthField
%       name of field that contains cell length (length_fitNew is current
%       one)

colorLetter     = upper(fluorcolor);
% rate
dcolorValueField = ['d' colorLetter '5'];
dcolorTimeField5 = ['d' colorLetter '5_time'];
% amount
colorSumField   = [colorLetter '5_sum'];
colorTimeField  = [colorLetter '_time'];
%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).(dcolorValueField)         = [];
    schnitzcells(i).(dcolorTimeField5)    = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    
    fluo_sum = s.(colorSumField);
    fluo_time = s.(colorTimeField);
    clear sD sE sP;
    
    %add parent. 
    %Hypothesis : at cell cutting, fluo is exactly shared between daughters
    % addition MW: according to ratio cell lengths. (TODO: maybe area would
    % be better even?)
    if s.P ~= 0
        sP = schnitzcells(s.P);
        if ~isempty(sP.(colorTimeField))
            % get ratio of length daughter:parent.
            thisGuysLengths = s.(lengthField);
            parentLengths = sP.(lengthField);
            lengthRatio = thisGuysLengths(1)/parentLengths(end);
            % add at start.
            fluo_sum = [sP.(colorSumField)(end)*lengthRatio fluo_sum];
                %fluo_sum = [sP.(colorSumField)(end)/2 fluo_sum]; % PN method
            fluo_time = [sP.(colorTimeField)(end) fluo_time];
        end
    end
    
    %add children
    if s.D ~= 0 && s.E ~= 0
        sD = schnitzcells(s.D);
        sE = schnitzcells(s.E);
        if ~isempty(sD.(colorTimeField)) && ~isempty(sE.(colorTimeField))
            % simply sum children's fluor 
            fluo_sum(end+1) = sD.(colorSumField)(1) + sE.(colorSumField)(1);
            fluo_time(end+1) = sD.(colorTimeField)(1);
        end
    end
    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           fit_window = [j-1 j j+1];
           xfit = fluo_time(fit_window);
           yfit = fluo_sum(fit_window);
           p = polyfit(xfit,yfit,1);
           
           s.(dcolorValueField)(end+1) = p(1);
           s.(dcolorTimeField5)(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


