function ns= neighbouringcells(img, cellno, dist)
% function ns= neighbouringcells(img, cellno, dist)
%
% returns neighbouring cells, lying less than dist from any edge of 
%  cell number cellno in img

if nargin == 2
    dist= 10;
end;

cell= +(img == cellno);
% extract subcell
[fx, fy]= find(cell);
extra= 5;
xmin= max(min(fx) - extra, 1);
xmax= min(max(fx) + extra, size(cell,1));
ymin= max(min(fy) - extra, 1);
ymax= min(max(fy) + extra, size(cell,2));
subcell= cell(xmin:xmax, ymin:ymax);

% dilate and find perimeter
subcell= imdilate(subcell, strel('disk', dist));
ns= [];
if length(find(cell>0))
    perim= bwperim(subcell);
    [px, py]= find(perim);
    px= px + xmin - 1;
    py= py + ymin - 1;

    % find unique neighbours
    for i= 1:length(px)
        ns= [ns img(px(i), py(i))];
    end;
    ns= unique(ns);
    ns(find(ns == 0))= [];
end