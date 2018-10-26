function [schnitzcells] = MW_fluorRate_anycolor(schnitzcells, fluorcolor, lengthField)
% Adds rate for any color
% DJK_addToSchnitzes_fluorRate_phase.
% 
% Later addition:
% Calculating the production rate has a few different approaches. Here we
% take d(N2-N1)/(area*dt) as the best approach (where N is the absolute
% signal, and dt the time difference between two signals N1 and N2). Note
% that this is different from d(C2-C1)/dt, which would take growth into
% account implicitly.
% Thus, 
% dC/dt = p - C * mu; where mu is growth and C concentration. 
% p can be determined by p = d(N2-N1)/(area*dt).
% This way of doing it is stored in the field:
% dcolorValueFieldDivAreaPx = ['d' colorLetter '5_divAreaPx'];
% 
% 
% input
%   - color    
%       letter that represents fluors color
%   - lengthField
%       name of field that contains cell length (length_fitNew is current
%       one)

colorLetter     = upper(fluorcolor);
% rate (where to put the rate)
dcolorValueField = ['d' colorLetter '5'];
dcolorValueFieldDivAreaPx = ['d' colorLetter '5_divAreaPx'];
dcolorTimeField5 = ['d' colorLetter '5_time'];
% amount (used to calculate the rate)
colorSumField   = [colorLetter '5_sum'];
colorTimeField  = [colorLetter '_time'];
%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).(dcolorValueField)    = [];
    schnitzcells(i).(dcolorValueFieldDivAreaPx)    = [];
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
    numelFrames = numel(fluo_sum);
    clear sD sE sP;
    
    % make a list of cell areas when fluor image was taken
    indicesFluor = ismember(s.time, fluo_time);
    areasAtFluor = s.areaPixels(indicesFluor);
    
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
            
            % dummy field for area to be consistent w. indices later
            areasAtFluor = [NaN areasAtFluor];
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
            % dummy field for area to be consistent w. indices later
            areasAtFluor = [areasAtFluor NaN];
        end
    end
    

    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           
           fit_window = [j-1 j j+1];
           xTofit = fluo_time(fit_window);
           yTofit = fluo_sum(fit_window);           
           p = polyfit(xTofit,yTofit,1);
           
           s.(dcolorValueField)(end+1) = p(1);
           s.(dcolorValueFieldDivAreaPx)(end+1) = p(1)/areasAtFluor(j);
           s.(dcolorTimeField5)(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


