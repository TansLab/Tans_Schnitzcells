function tree = AscendField(schni,fieldName,n)
    s = schni(n);
    sP = s.P;
     
    if sP == 0 
        tree = [];
    else
        sparent = schni(sP);
        tree = [AscendField(schni,fieldName,sP) ; sparent.(fieldName)(:)];
    end
    
end