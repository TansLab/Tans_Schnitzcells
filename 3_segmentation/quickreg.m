function [xsubreg, xshift, xback, xbinning] = quickreg(L, imx, rect, deltamax, fullsize)
% function [xsubreg, xshift, xback, xbinning] = quickreg(L, imx, rect, deltamax, fullsize)
% 
%  calculates translation between phase and fluorescent images

% load images if necessary

if min(size(imx))==1,
    if isempty(findstr('.tif', imx)),
        imx = [imx,'.tif'];
    end;
    if exist(imx)==2
        imx= imread(imx);
    end
end

xbinning = 1;
if size(imx)==round(0.5*fullsize),
    disp('your fluor image is 2x2 binned...');
    imx = imresize(imx,2,'nearest');
    xbinning = 2;
end;

if min(size(imx))>1
% get subimages
imx1= double(imx(rect(1):rect(3), rect(2):rect(4)));
LL= +(L(1+deltamax:end-deltamax, 1+deltamax:end-deltamax) > 0);
% try all possible translations
for di= -deltamax:deltamax,
    for dj= -deltamax:deltamax,
        imx2= imx1(deltamax+di+1:end-deltamax+di, deltamax+dj+1:end-deltamax+dj);
        xsum(di+deltamax+1,dj+deltamax+1)= sum(imx2(LL > 0));
    end;
end;
% best translation is the one with largest csum (the most white pixels in LL)
[xbesti, xbestj]= find(xsum == max2(xsum));
xbesti= xbesti - deltamax - 1;
xbestj= xbestj - deltamax - 1;
% record translation
xshift= [xbesti(1) xbestj(1)];
% record translated images
xsubreg= imx((rect(1):rect(3))+xbesti(1), (rect(2):rect(4))+xbestj(1));
% calculate background fluorescence
imxb = double(imx);
imxb(rect(1):rect(3), rect(2):rect(4))=0;
imxbvect=imxb(imxb>0);
xback=median(imxbvect);
else
    xsubreg=[];
    xshift=[];
    xback=[];
end;
