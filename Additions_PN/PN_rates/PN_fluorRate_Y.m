function [schnitzcells] = PN_fluorRate_Y(schnitzcells)

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).dY5         = [];
    schnitzcells(i).dY5_time    = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    
    fluo_sum = s.Y5_sum;
    fluo_time = s.Y_time;
    clear sD sE sP;
    
    %add parent. 
    %Hypothesis : at cell cutting, fluo is exactly shared between daughters
    if s.P ~= 0
        sP = schnitzcells(s.P);
        if ~isempty(sP.Y_time)
            fluo_sum = [sP.Y5_sum(end)/2 fluo_sum];
            fluo_time = [sP.Y_time(end) fluo_time];
        end
    end
    
    %add children
    if s.D ~= 0 && s.E ~= 0
        sD = schnitzcells(s.D);
        sE = schnitzcells(s.E);
        if ~isempty(sD.Y_time) && ~isempty(sE.Y_time)
            fluo_sum(end+1) = sD.Y5_sum(1) + sE.Y5_sum(1);
            fluo_time(end+1) = sD.Y_time(1);
        end
    end
    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           fit_window = [j-1 j j+1];
           xfit = fluo_time(fit_window);
           yfit = fluo_sum(fit_window);
           p = polyfit(xfit,yfit,1);
           
           s.dY5(end+1) = p(1);
           s.dY5_time(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


