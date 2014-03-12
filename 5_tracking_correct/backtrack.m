function [cpred, sc]= backtrack(today, yesterday, cellno, trans, transmax, len, padsize)
% function [cpred, sc]= backtrack(today, yesterday, cellno, trans, transmax, len, padsize)
%
% backward tracking algorithm; tracks cellno in todays image to the equivalent cell in 
%  yesterdays image
% Called from trackcomplete.m

spurpts= 5;  % number of spurs removed from thin to remove forked ends
sc= 0;

% find centroids
cell= +(today == cellno);
stats= regionprops(cell, 'Centroid');
cen= round([stats.Centroid]);
if length(cen)<2
    cpred= 1; %was =-1;
    return;
else
x= cen(2);
y= cen(1);
% check centroid lies inside cell
if cell(x, y) ~= 1
    disp(' centroid outside cell');
    % extract subcell
    [fx, fy]= find(cell);
    extra= 5;
    xmin= max(min(fx) - extra, 1);
    xmax= min(max(fx) + extra, size(cell,1));
    ymin= max(min(fy) - extra, 1);
    ymax= min(max(fy) + extra, size(cell,2));
    subcell= cell(xmin:xmax, ymin:ymax);
    
    % find thin
    thin= bwmorphmelow(subcell, 'thin', inf);
    thins= bwmorphmelow(thin, 'spur', spurpts);
    [px, py]= walkthin(thins);
    
    % centroid is center point on thin
    cen= round(length(px)/2);
    x= px(cen) + xmin - 1;
    y= py(cen) + ymin - 1;
end;


% extract subimages
[todaypadnew,yesterdaypadnew]=padandshift(today,yesterday,[0,0],len/2-1,0);
tsub= todaypadnew(x:x + len, y:y + len);
ysub = yesterdaypadnew(padsize + x + trans(1):padsize + x + trans(1) + len,...
    padsize + y + trans(2):padsize + y + trans(2) + len);
% tsub= today(x - len/2:x + len/2, y - len/2:y + len/2);
% ysub = yesterday(padsize + x + trans(1) - len/2:padsize + x + trans(1) + len/2,...
%     padsize + y + trans(2) - len/2:padsize + y + trans(2) + len/2);

% calculate offset
t= imoffset(tsub > 0, ysub > 0);
if max(abs(t)) < min([len/2 transmax]),
    xn= len/2 + t(1);
    yn= len/2 + t(2);
else    
    disp([' ignoring offset; tmax= ',num2str(max(abs(t)))]);
    sc= 1;
    xn= len/2;
    yn= len/2;
end;


% extract cell number from yesterday's image
cpred= ysub(xn, yn);
if cpred==0 
    % translated centroid lies outside cells
    sc= 1;
    
    % find labels of cells
    u= unique(ysub(:));
    u(1)= [];
    
    % find distance of each cell from (xn_ye, yn_ye)
    d= [];
    if isempty(u)
        cpred= 1; %was =-1;
    else
        for k= 1:length(u)
            d= [d mind(xn, yn, ysub, u(k))];
        end;
        % find closest cell to (xn_ye, yn_ye)
        cpred= u(minI(d));
    end;
end;
end