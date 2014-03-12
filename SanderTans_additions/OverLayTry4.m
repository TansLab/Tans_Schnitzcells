cmap=repmat([0 0 0;0 0 0],[200,1]);
%cmap=repmat([0 0 0;0.4 0.4 0.4],[200,1]);
miny=min(min(yreg));
maxy=max(max(yreg));
nrvalues=maxy-miny;
yregN = yreg-miny;
yregRGB=ind2rgb(yregN,jet(double(nrvalues)));

SE=strel('square',3);

perim = LNsub-imerode(LNsub,SE);
perimRGB = ind2rgb(perim,cmap);
%figure;image(perimRGB);

immask = ~(perim>0);

immaskRGB=repmat(double(immask),[1,1,3]);
%figure;image(immaskRGB);
combRGB = (yregRGB .*  (immaskRGB)) + perimRGB ;
figure;
image(combRGB);
