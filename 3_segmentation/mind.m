
function dist = mind(pt_x, pt_y, img, cellno);
% function dist = mind(pt_x, pt_y, img, cellno);
%
% returns minimum distance of (pt_x, pt_y) from cell denoted by cellno
%  in image img

[fy,fx] = find(img==cellno);
dy = fy - pt_y;
dx = fx - pt_x;
d = sqrt(dx.^2 + dy.^2);
dist = min(d);

