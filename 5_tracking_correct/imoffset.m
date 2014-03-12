function trans= imoffset(im1,im2);
% function trans= imoffset(im1,im2);
%
% takes two images calculates translation between them
% using cross correlation. 
% Images must be of same size.

% calculate normalised cross correlation
corr= normxcorr2(im1,im2);
[maxc,i]= max(abs(corr(:)));
[ypeak, xpeak]= ind2sub(size(corr),i(1));

% calculate offset 
trans= [(xpeak-size(im1,2)) (ypeak-size(im1,1))];
