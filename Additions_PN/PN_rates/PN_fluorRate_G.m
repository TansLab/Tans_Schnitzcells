function [schnitzcells] = PN_fluorRate_G(schnitzcells)
% same as PN_fluorRate but for GFP instead of YFP 

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).dG5         = [];
    schnitzcells(i).dG5_time    = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    
    fluo_sum = s.G5_sum;
    fluo_time = s.G_time;
    clear sD sE sP;
    
    %add parent. 
    %Hypothesis : at cell cutting, fluo is exactly shared between daughters
    if s.P ~= 0
        sP = schnitzcells(s.P);
        if ~isempty(sP.G_time)
            fluo_sum = [sP.G5_sum(end)/2 fluo_sum];
            fluo_time = [sP.G_time(end) fluo_time];
        end
    end
    
    %add children
    if s.D ~= 0 && s.E ~= 0
        sD = schnitzcells(s.D);
        sE = schnitzcells(s.E);
        if ~isempty(sD.G_time) && ~isempty(sE.G_time)
            fluo_sum(end+1) = sD.G5_sum(1) + sE.G5_sum(1);
            fluo_time(end+1) = sD.G_time(1);
        end
    end
    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           fit_window = [j-1 j j+1];
           xfit = fluo_time(fit_window);
           yfit = fluo_sum(fit_window);
           p = polyfit(xfit,yfit,1);
           
           s.dG5(end+1) = p(1);
           s.dG5_time(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


