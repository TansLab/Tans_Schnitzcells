function imnew = DJK_randomRenumberImage(imold);
% DJK: copied from renumberimage and added randomness

imnew = imold;
u = unique(imold(:));
r = randperm(length(u)-1);

% u(1) = 0 is background, so start with u(2)
for i = 2:length(u),
  imnew(imold==u(i))=r(i-1);
end;

if islogical(imnew),
  imnew = +imnew;
end;
