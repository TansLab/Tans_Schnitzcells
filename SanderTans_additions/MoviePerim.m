cmap=repmat([0 0 0;0 0 0],[200,1]); %first color is background (index=1), followed by cell nr index
%cmap=repmat([0 0 0;0.4 0.4 0.4],[200,1]);

figure;
image(yreg,'CDataMapping', 'scaled');


miny=min(min(yreg))
maxy=max(max(yreg))
nrvalues=maxy-miny
% yregN = yreg*20;
yregN = (yreg-miny+1)*1000;
yregRGB=ind2rgb(yregN,jet(double(nrvalues)));
% yregRGB=ind2rgb(yregN,jet);

% figure;
% image(yregRGB,'CDataMapping','scaled');

SE=strel('square',3);

perim = Lc-imerode(Lc,SE);
%perim = LNsub-imerode(LNsub,SE);
perimRGB = ind2rgb(perim,cmap);
%figure;image(perimRGB);

immask = ~(perim>0);

immaskRGB=repmat(double(immask),[1,1,3]);
%figure;image(immaskRGB);
combRGB = (yregRGB .*  (immaskRGB)) + perimRGB ;
figure;
image(combRGB);
