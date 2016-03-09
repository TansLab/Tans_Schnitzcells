
function [px, py] = walkthin(thin)
% function [px, py] = walkthin(thin)
%
% returns ordered list of points along thin
%   starting from a spur point

% add empty border
thin2= zeros(size(thin)+2);
thin2(2:end-1, 2:end-1)= thin;
nopts= length(find(thin > 0));

if nopts < 5
    px= [];
    py= [];
    return
end;

% find spurs
spur= thin2 & ~bwmorphmelow(thin2,'spur',1);
[sx, sy]= find(spur);
if isempty(sx)
    px= [];
    py= [];
    return
end;

% find path
i= 1; 
px(1)= sx(1);
py(1)= sy(1);
while any2(thin2),
    
    thin2(px(i), py(i))= 0;
    [nx, ny]= neighbours(thin2, px(i), py(i));
    
    i= i+1;
    if length(nx) > 1
        % pick closest
        dI= minI(sqrt((nx - px(i-1)).^2 + (ny - py(i-1)).^2));
        px(i)= nx(dI);
        py(i)= ny(dI);
    elseif length(nx) == 1
        px(i)= nx(1);
        py(i)= ny(1);
    else
        break;
    end;
    
    if i > nopts
        break;
    end;
        
end;

px= px - 1;
py= py - 1;
    