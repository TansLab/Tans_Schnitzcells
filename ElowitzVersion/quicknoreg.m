function [xsubreg, xshift, xback, xbinning] = quicknoreg(L, imx, rect, deltamax, fullsize)
% function [xsubreg, xshift, xback, xbinning] = quickreg(L, imx, rect, deltamax)
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


sizeratio = fullsize./size(imx);
if sizeratio(1)==sizeratio(2) %& sizeratio(1)==round(sizeratio(1))
    if sizeratio(1)~=1
        disp(['fluor image is ',...
            num2str(sizeratio(1)),'x',num2str(sizeratio(1)),' binned.']);
        imx = imresize(imx,sizeratio(1),'nearest');
    end
    xbinning = sizeratio(1);
else
    disp('fluor image dimensions have different proportions that phase image.');
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
