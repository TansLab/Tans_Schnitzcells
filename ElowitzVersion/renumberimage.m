
function imnew = renumberimage(imold);

imnew = imold;
u = unique(imold(:));
for i = 2:length(u),
    imnew(imold==u(i))=i-1;
end;

if islogical(imnew),
    imnew = +imnew;
end;
