function bw2 = bwmorphmelow(bw1,operation,n);

marg = 3;

[fx,fy] = find(bw1);
x1 = max(min(fx)-marg,1);
x2 = min(max(fx)+marg,size(bw1,1));

y1 = max(min(fy)-marg,1);
y2 = min(max(fy)+marg,size(bw1,2));

if nargin <= 2,
    Z = bwmorph(bw1(x1:x2,y1:y2),operation);
else
    Z = bwmorph(bw1(x1:x2,y1:y2),operation,n);
end;

bw2 = zeros(size(bw1));
bw2(x1:x2,y1:y2) = Z;

