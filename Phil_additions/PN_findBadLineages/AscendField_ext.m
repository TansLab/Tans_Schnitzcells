function tree = AscendField(schni,fieldName,n)

% NOT MODIFIED YET!!! DONT USE !!!
    s = schni(n);
    sP = s.P;
     
    if sP == 0 
        tree = [];
    else
        sparent = schni(sP);
        tree = [AscendField(schni,fieldName,sP) ; sparent.(fieldName)(:)];
    end
    
end