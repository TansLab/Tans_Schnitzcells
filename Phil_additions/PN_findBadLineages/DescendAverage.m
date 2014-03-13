function avTree = DescendAverage(schni,fieldName,n)
    s = schni(n);

    if s.D == 0 || s.E == 0
        avTree = [];
    else
        
        kidD = schni(s.D);
        kidE = schni(s.E);
        kD = [kidD.(fieldName)(:) ; DescendAverage(schni,fieldName,s.D)];
        kE = [kidE.(fieldName)(:) ; DescendAverage(schni,fieldName,s.E)];
        avTree = mean([kD kE],2);

    end
    
end