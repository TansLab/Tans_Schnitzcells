
function [nx, ny, colour]= neighbours(img, x, y, bsize)
% function [nx, ny, colour]= neighbours(img, x, y, bsize)
%
% returns neighbours (and their colours) of pt(x,y) in img
% i.e. all points lying in the edge of a square of length 2*bsize

if nargin == 3
    bsize= 1;
end;
colour= [];

% pad image
img2= zeros(size(img) + 2*bsize);
img2(bsize+1:size(img2,1)-bsize, bsize+1:size(img2,2)-bsize)= img;

% extract subimage
bxmin= x;
bxmax= x + 2*bsize;
bymin= y;
bymax= y + 2*bsize;
subimg= img2(bxmin:bxmax, bymin:bymax);
if bsize > 1
    subimg(2:2*bsize, 2:2*bsize)= zeros(2*bsize-1);
end;
subimg(bsize+1, bsize+1)= 0;

% find neighbours
[nx, ny]= find(subimg > 0);
% find pixel values of neighbours
for i= 1:length(nx)
    colour(i,1)= subimg(nx(i), ny(i));
end;

% translate back to original image
nx= nx + bxmin - 1 - bsize;
ny= ny + bymin - 1 - bsize;
