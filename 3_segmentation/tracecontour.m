

function [px,py,con,conim] = tracecontour(im0);

% expand the image by one blank pixel.

N = [ 1 2 3 ; 8 0 4 ; 7 6 5];

NN(1,1:2) = [-1 -1];
NN(2,1:2) = [0  -1];
NN(3,1:2) = [1  -1];
NN(4,1:2) = [1  0];
NN(5,1:2) = [1  1];
NN(6,1:2) = [0 1];
NN(7,1:2) = [-1 1];
NN(8,1:2) = [-1 0];

im = zeros(size(im0)+2);
im(2:end-1,2:end-1) = im0;

con = [];
[fx,fy] = find(im);          % find initial point.
if length(fx)==0
    px=[];
    py=[];
    conim=[];
    return;
end
p1 = [fx(1) fy(1)];          % set current point p1
con(1,1:2) = [p1(1) p1(2)];  % add initial point to list.
p0 = p1-1;                   % initial departure point, p0
s = p1;

% algorithm: start at current point, p1. 
% try all the points around it, in a clockwise direction, starting with the
% point that comes after p0.
% eventually you will find an "on" point.
% make this point your new p1 and make the previous p1 your new p0.
% repeat.
% stop when the initial p1 point is in your exploration neighberhood.

d = p0-p1;
p0num = N(2+d(1),2+d(2));
dontstop = 0;
while dontstop<1,
    
    nextp0num = 1+mod(p0num,8);
    [f1,f2] = find(N==nextp0num);
    nextp0 = p1 + [f1(1)-2 f2(1)-2];
    
    if im(nextp0(1),nextp0(2)),
        p0 = p1;
        p1 = nextp0;
        d = p0-p1;
        con(end+1,1:2) = p1;
        p0num = N(2+d(1),2+d(2));
    else
        p0num = nextp0num;
    end;
    if nextp0 == s;
        dontstop = dontstop +  1;
    end;
end;

con = con - 1;
if nargout>1,
    conim = zeros(size(im0));
    for i = 1:size(con,1),
        conim(con(i,1),con(i,2)) = 1;
    end;
end;

       
px = con(:,1);
py = con(:,2);

