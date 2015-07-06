function [schnitzcells] = PN_fluorRate_C(schnitzcells)
% same as PN_fluoRate but for CFP instead of YFP 
% MW TODO: generalize these functions same as
% DJK_addToSchnitzes_fluorRate_phase.

%--------------------------------------------------------------------------
% INITIALIZE NEW MEASUREMENTS
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    schnitzcells(i).dC5         = [];
    schnitzcells(i).dC5_time    = [];
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
    s = schnitzcells(i);
    
    fluo_sum = s.C5_sum;
    fluo_time = s.C_time;
    clear sD sE sP;
    
    %add parent. 
    %Hypothesis : at cell cutting, fluo is exactly shared between daughters
    if s.P ~= 0
        sP = schnitzcells(s.P);
        if ~isempty(sP.C_time)
            fluo_sum = [sP.C5_sum(end)/2 fluo_sum];
                % MW '15/04 TODO possible bug: I think this should be
                % scaled w. the length of the daughter cells.
            fluo_time = [sP.C_time(end) fluo_time];
        end
    end
    
    %add children
    if s.D ~= 0 && s.E ~= 0
        sD = schnitzcells(s.D);
        sE = schnitzcells(s.E);
        if ~isempty(sD.C_time) && ~isempty(sE.C_time)
            fluo_sum(end+1) = sD.C5_sum(1) + sE.C5_sum(1);
            fluo_time(end+1) = sD.C_time(1);
        end
    end
    
    %add a rate computed linearly fitted over 3 points and centered at the middle one
    if length(fluo_sum)>2
       for j = 2:length(fluo_sum)-1
           fit_window = [j-1 j j+1];
           xfit = fluo_time(fit_window);
           yfit = fluo_sum(fit_window);
           p = polyfit(xfit,yfit,1);
           
           s.dC5(end+1) = p(1);
           s.dC5_time(end+1) = fluo_time(j);
       end
    end

    schnitzcells(i) = s;
end


