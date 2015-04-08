function [xsubreg, xshift, xback, xbinning] = quicknoreg(L, imx, rect, deltamax, fullsize)
% function [xsubreg, xshift, xback, xbinning] = quickreg(L, imx, rect, deltamax)
%
% - Calculates translation between phase and fluorescent images
% - Also calculates binning use for fluor image based on fluor image ("imx" 
% gives the fluor file and is loaded), and the given size of the phase
% image ("fullsize").
%
% load images if necessary

if min(size(imx))==1,
    if isempty(findstr('.tif', imx)),
        imx = [imx,'.tif'];
    end;
    if exist(imx)==2
        imx= imread(imx);
    end
end


sizeratio = fullsize./size(imx);
% I've build the following warning in because otherwise the error
% introduced by this issue would remain hidden. There might be prettier
% ways to do this though. 2048 is the current resolution of the camera,
% when taking pictures at binning=1. MW 2015/04
CURRENTCAMERARESOLUTION = 2048;
if fullsize(1)<CURRENTCAMERARESOLUTION 
    disp(['WARNING: size of your phase image is smaller than what is typical (per 2015);'...
          ' If you have used binning for your phase image, binning of the fluor image will be determined incorrectly.']);
end
if sizeratio(1)==sizeratio(2) %& sizeratio(1)==round(sizeratio(1))
    if sizeratio(1)~=1
        disp(['fluor image is ',...
            num2str(sizeratio(1)),'x',num2str(sizeratio(1)),' binned.']);
        imx = imresize_old(imx,sizeratio(1),'nearest');
    end
    xbinning = sizeratio(1);
else
    disp('fluor image dimensions have different proportions than phase image.');
    error(' not equipped for such cases.');
end;
% end nitzan editing 2005June25

if min(size(imx))>1
    % get subimages
    imx1= double(imx(rect(1):rect(3), rect(2):rect(4)));
    LL= +(L(1+deltamax:end-deltamax, 1+deltamax:end-deltamax) > 0);
    % best translation is the one with largest csum (the most white pixels in LL)

    xshift= [0 0];
    % record translated images
    %xsubreg= imx((rect(1):rect(3))+xbesti(1), (rect(2):rect(4))+xbestj(1));
    % JCR: Fixed bounds problems
    xsubreg = imx(max(rect(1),1):min(rect(3),size(imx,1)),max(1,rect(2)):min(size(imx,2),rect(4)));
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
