function suspicous = MW_getsuspicous(s)

    suspicous=[];

    for idx = 1:length(s)
        if any(s(idx).muP11_all < 0)
            suspicous = [suspicous, idx];
        end
    end

end