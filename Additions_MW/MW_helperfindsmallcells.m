

allCellNrs=unique(Lc(:))';
for cellNr = allCellNrs
    
    if sum(Lc(:)==cellNr) < 100
       cellNr 
    end
    
end