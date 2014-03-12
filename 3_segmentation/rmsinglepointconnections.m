function cell2= rmsinglepointconnections(cell)
% function cell2= rmsinglepointconnections(cell)
%
% removes single pixels that connect otherwise disparate regions of image

cell= +(cell > 0);
thin= bwmorphmelow(cell, 'thin', inf);
perim= bwperim(cell);
% potential points to be removed
[fx, fy]= find(thin & perim);

cell2= cell;
for i= 1:length(fx)
    
    % extract subimage
    extra= 5;
    xmin= max(min(fx(i)) - extra, 1);
    xmax= min(max(fx(i)) + extra, size(cell,1));
    ymin= max(min(fy(i)) - extra, 1);
    ymax= min(max(fy(i)) + extra, size(cell,2));
    subcell= cell(xmin:xmax, ymin:ymax);
    
    [ell, nbef]= bwlabel(subcell, 4);
    subcell(extra+1, extra+1)= 0;
    [ell, naft]= bwlabel(subcell, 4);
    if naft > nbef
        cell2(fx(i), fy(i))= 0;
    end;
    
end;

cell2= bwlabel(cell2, 4);

