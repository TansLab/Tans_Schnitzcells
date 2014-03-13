function listkids = Descend(schni,n)
%includes the initial schnitz

s = schni(n);
    sD = s.D;
    sE = s.E;
    
    if s.D == 0 && s.E == 0
        listkids = [n];
    elseif s.D == 0
        listkids = [n Descend(schni,sE)];
    elseif s.E == 0
        listkids = [n Descend(schni,sD)];
    else
        listkids = [n Descend(schni,sD) Descend(schni,sE)];
    end
    
end