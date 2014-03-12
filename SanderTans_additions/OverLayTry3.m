
[zx,zy]=size(yreg);
blank=zeros(zx,zy);
blankRGB=ind2rgb(blank,gray(2));


figure;
colormap('jet');
image(yregN,'CDataMapping','scaled');
caxis([0 nrvalues]);

miny=min(min(yreg));
maxy=max(max(yreg));
nrvalues=maxy-miny;
yregN = yreg-miny;
yregRGB=ind2rgb(yregN,jet(double(nrvalues)));
figure;image(yregRGB);
SE=strel('square',3);
immask=~(LNsub-imerode(LNsub,SE)>0);
immaskRGB=repmat(double(immask),[1,1,3]);
combRGB = (yregRGB .*  (immaskRGB)) + ~immaskRGB ;
figure;
image(combRGB);

figure;
image(immask,'CDataMapping','scaled');
colormap('gray');

figure;
immaskRGB=ind2rgb(immask,gray(2));
image(immaskRGB);

cmap = colormap(jet(double(nrvalues)));
yregRGB=ind2rgb(yregN,cmap);

figure;
imy=image(yregRGB);
hold on;
imm=image(blankRGB);
set(imm,'AlphaData',double(immask));
result = getimage(imgca);

yregRGB(immask>0)=[1 1 1];
repl=( bRGB((1:size(aRGB,1)),(1:size(aRGB,2)),:) .*  double(~aRGB) )+(aRGB .* double(aRGB));
repl= bRGB .*  (~aRGB) + aRGB ;
repl= bRGB .*  (~aRGB) + cRGB ;

colormap('gray');
hold on;
imm=image(90*immask);
colormap('gray');
