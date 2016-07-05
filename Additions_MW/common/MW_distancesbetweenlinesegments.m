function distancesxy = MW_distancesbetweenlinesegments(x,y)

% euclidian distance is sqrt((x2-x1)^2+(y2-y1)^2)

dx = x(2:end)-x(1:end-1);
dy = y(2:end)-y(1:end-1);

distancesxy = sqrt(dx.^2+dy.^2);

end