function [schnitzcells] = PN_fluorRate_R(schnitzcells)
% same as PN_fluorRate but for RFP instead of YFP 

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).dR5         = [];
    schnitzcells(i).dR5_time    = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    
    fluo_sum = s.R5_sum;
    fluo_time = s.R_time;
    clear sD sE sP;
    
    %add parent. 
    %Hypothesis : at cell cutting, fluo is exactly shared between daughters
    if s.P ~= 0
        sP = schnitzcells(s.P);
        if ~isempty(sP.R_time)
            fluo_sum = [sP.R5_sum(end)/2 fluo_sum];
            fluo_time = [sP.R_time(end) fluo_time];
        end
    end
    
    %add children
    if s.D ~= 0 && s.E ~= 0
        sD = schnitzcells(s.D);
        sE = schnitzcells(s.E);
        if ~isempty(sD.R_time) && ~isempty(sE.R_time)
            fluo_sum(end+1) = sD.R5_sum(1) + sE.R5_sum(1);
            fluo_time(end+1) = sD.R_time(1);
        end
    end
    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           fit_window = [j-1 j j+1];
           xfit = fluo_time(fit_window);
           yfit = fluo_sum(fit_window);
           p = polyfit(xfit,yfit,1);
           
           s.dR5(end+1) = p(1);
           s.dR5_time(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


