function [newcell, sc]= forwardtrack(today, yesterday, cellno, trans, transmax, len, padsize)
% function [newcell, sc]= forwardtrack(today, yesterday, cellno, trans, transmax, len, padsize)
%
% forward tracking algorithm; tracks cellno in yesterdays image to the equivalent cell or cells in 
%  todays image. Uses two `centroids', 1/4 and 3/4 along thin
% Called from trackcomplete.m



spurpts= 5;  % number of spurs removed from thin to remove forked ends
sc= 0;

% find centroids
cell= +(yesterday == cellno);
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
if length(find(thins > 0)) < 5
    [px, py]= walkthin(thin);
else
    [px, py]= walkthin(thins);
end;

if length(px)>0 % added by nitzan 2002-11-21
% centroids is are thirds of thin
cen(1)= round(length(px)/4);
cen(2)= round(3*length(px)/4); 
for k= 1:2
    x(k)= px(cen(k)) + xmin - 1;
    y(k)= py(cen(k)) + ymin - 1;
end;

% check each centroid in turn
for k= 1:2
    
    % extract subimages
    [todaypadnew,yesterdaypadnew]=padandshift(today,yesterday,[0,0],len/2-1,0);
    ysub{k}= yesterday(x(k):x(k) + len, y(k):y(k) + len);
    tsub{k}= today(padsize + x(k) - trans(1):padsize + x(k) - trans(1) + len,...
        padsize + y(k) - trans(2):padsize + y(k) - trans(2) + len);
%     ysub{k}= yesterday(x(k) - len/2:x(k) + len/2, y(k) - len/2:y(k) + len/2);
%     tsub{k}= today(padsize + x(k) - trans(1) - len/2:padsize + x(k) - trans(1) + len/2,...
%         padsize + y(k) - trans(2) - len/2:padsize + y(k) - trans(2) + len/2);
    
    % calculate offset
    t= imoffset(ysub{k} > 0, tsub{k} > 0);
    if max(abs(t)) < min([len/2 transmax])
        xn{k}= len/2 + t(1);
        yn{k}= len/2 + t(2);
    else    
        disp([' ignoring offset; tmax= ',num2str(max(abs(t)))]);
        sc= 1;
        xn{k}= len/2;
        yn{k}= len/2;
    end;
    
    
    % extract cell number from yesterday's image
    cpred{k}= tsub{k}(xn{k}, yn{k});
    if cpred{k}==0  
        % translated centroid lies outside cells
        sc= 1;
                
        % find labels of cells
        u= unique(tsub{k}(:));
        u(1)= [];
        
        % find distance of each cell from (xn_ye, yn_ye)
        d= [];
        if isempty(u)
            cpred{k}= -1;
        else
            for j= 1:length(u)
                d= [d mind(xn{k}, yn{k}, tsub{k}, u(j))];
            end;
            % find closest cell to (xn_ye, yn_ye)
            cpred{k}= u(minI(d));
        end;
    end;
    
end;


% amalgamate results
if cpred{1} == cpred{2}
    newcell= cpred{1};
else
    newcell= [cpred{1} cpred{2}];
end;
% the following 4 lines were added by nitzan 2002-11-21:
else
newcell=0;
sc=-1;
end    