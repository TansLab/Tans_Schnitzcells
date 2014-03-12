function dL = carefuldilate(L,blob,niter,refimage);
% function dL = carefuldilate(L, blob, niter, refimage)
%
% dilates L by strel blob niter times
% if regfimage is included, points are only added that exist in refimage

if nargin < 4,
    refimage = ones(size(L));
end;
    

for i = 1:niter,
    
    LD = imdilate(L,blob);
    BW0 = (L>0);
    BW1 = imdilate(BW0,blob);
    NewPix = BW1 & ~BW0 & refimage;
    
    L(NewPix) = LD(NewPix); 
    
end;
  
dL = L;
