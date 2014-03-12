function Iout = PN_reseed(Iin,phim,ppx,ppy)

Iout = Iin;
edgeIm = edge(phim,'log',0,2); 
newCell = imfill(edgeIm,[ppx ppy]) & ~edgeIm;

intersectionIm = newCell & Iout;

if max2(intersectionIm) > 0 
    localValue = max2(Iout(logical(intersectionIm)));
    Iout(newCell) = localValue; 
else
    Iout(newCell) = max2(Iin)+1;
end


end