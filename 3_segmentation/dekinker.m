% 
function cutcell= dekinker(cell, radius, mincelllength, angthresh)
% function cutcell= dekinker(cell, radius, mincelllength, angthresh)
%
% cuts kinky cells by finding the angle between two points on thin a distance
%   radius away from a centre point
% cuts points whose angle is a local minimum along thin
% angthresh:        larger this is the more points are cut
% mincelllength:    cut will be ignored if it creates 'cells' smaller than this

angle=[]; % added by nitzan 21/8/2002, to take care of annoying warning
          % of uninitialized variable "angle" at line 92.
          
mindist= 5; % distance from pt used to test if it is a minimum
diffthresh= 3; % size of jump in diff of minpts  

% make sure image is black and white and not logical
cellno= max2(cell);
cell= +(cell > 0);
% cell= rmsinglepointconnections(cell);

% extract subimages
[fx, fy]= find(cell);
extra= 5;
xmin= max(min(fx) - extra, 1);
xmax= min(max(fx) + extra, size(cell,1));
ymin= max(min(fy) - extra, 1);
ymax= min(max(fy) + extra, size(cell,2));
subcell= cell(xmin:xmax, ymin:ymax);
originalsubcell= subcell;

% find thin
thin0 = bwmorphmelow(subcell, 'thin', inf);
% reduce so that there is only two spurs
thin= thin0;
spur= thin & ~bwmorphmelow(thin, 'spur', 1);
while length(find(spur)) > 2
    thin= thin & ~spur;
    spur= thin & ~bwmorphmelow(thin, 'spur', 1);
end;



% FIND ANGLES ALONG THIN
[fx,fy]= walkthin(thin);
for i= 1:length(fx)
    
    % find distances
    dist= sqrt((fx(i)-fx).^2 + (fy(i)-fy).^2);
    
    % find thetas (between 0 and 2 pi)
    thetas= pi + atan2(fy(i)-fy, fx(i)-fx);
    for j= 1:length(thetas),
        if thetas(j) < 0,
            thetas(j)= thetas(j) + 2*pi;
        end;
    end;
    
    % find points on thin between two distances
    pts{i} = find((dist > 0.75*radius) & (dist < 1.4*radius));
    if isempty(pts{i})
        
        disp(['Dekinker: no pixels near radius. Cell number: ', num2str(cellno)]);
        cutcell= bwlabel(cell, 4);
        return;
        
    else
        
        thets{i}= thetas(pts{i});
        
        % calculate thets from some arbitrary thet
        arbthet= thets{i}(1);
        dthet= [];
        for j= 2:length(thets{i}),
            dthet(j)= abs(thets{i}(j) - arbthet);
            if dthet(j) > 2*pi,
                dthet(j)= dthet(j) - 2*pi;
            elseif dthet(j) < 0,
                dthet(j)= dthet(j) + 2*pi;
            elseif dthet(j) > pi,
                dthet(j)= 2*pi-dthet(j);
            end;
            % require max angle difference between 
            % points to left of thin point and points to right
            angle(i)= max(dthet);
        end;
        
    end;
    
end;


% FIND LOCAL MINIMA IN ANGLE (potential points to be cut)
minpts= [];
for i= 1 + mindist:length(angle) - mindist
    if (angle(i - mindist) > angle(i)) & (angle(i + mindist) > angle(i))
        minpts= [minpts i];
    end;
end; 
if ~isempty(minpts)
    % minpts contains multiple points for each minimum
    % find boundaries between sets of points
    minboundaries= unique([1 find(diff(minpts) > diffthresh) length(minpts)]);
    
    % each minimum is the average of each set of points associated with it
    if length(minboundaries)> 1
        
        for i= 2:length(minboundaries)
            subminpts= minpts(minboundaries(i-1) + 1:minboundaries(i));
            ami= minI(angle(subminpts));
            mins(i-1)= subminpts(ami);
        end;
        
        % choose large minima only
        kinkpts= [mins(find(angle(mins) < angthresh))];
    else
        kinkpts= [];
    end;
    
else
    kinkpts= [];
end;



% COMBINE NEIGHBOURING POINTS
if isempty(kinkpts)
    kinkx= [];
    kinky= [];
else
    kinkx= fx(kinkpts);
    kinky= fy(kinkpts);
end;

i= 1;
while i < length(kinkx),
    
    d= sqrt((kinkx(i)-kinkx(i+1))^2 + (kinky(i)-kinky(i+1))^2);        
    if d < mincelllength
        kinkx(i)= round(mean([kinkx(i) kinkx(i+1)]));
        kinky(i)= round(mean([kinky(i) kinky(i+1)]));
        
        kinkx(i+1)= [];
        kinky(i+1)= [];
    else
        i= i + 1;
    end;
    
end;


% CUT CELLS
% now divide the cell into seperate cells by cutting it across
perim= bwperim(imdilate(subcell, strel('disk',1)));

for i= 1:length(kinkx)  
        
    % extract box around kink
    bsize= 8;
    sxmin= max(1, kinkx(i) - bsize);
    sxmax= min(size(perim,1), kinkx(i) + bsize);
    symin= max(1, kinky(i) - bsize);
    symax= min(size(perim,2), kinky(i) + bsize);
    subperim= perim(sxmin:sxmax, symin:symax);
    [subperim, noperims]= bwlabel(subperim);
    
    if any(any(subperim)) & (noperims ~= 1),
        
        % kink is not near end of cell. Go ahead and cut.
        currcell= subcell;
        [px, py]= find(subperim> 0);
        
        % find distances to perimeter from cutpt
        d= sqrt((px - bsize - 1).^2 + (py - bsize - 1).^2);
        [ds, di]= sort(d);
        
        % find first cutting point on perimeter
        cutperim1x= px(di(1)) + kinkx(i) - bsize - 1;
        cutperim1y= py(di(1)) + kinky(i) - bsize - 1;
        colour1= subperim(px(di(1)), py(di(1)));
        
        % find second cutting point on perimeter
        colour= colour1;
        j= 2;
        while colour == colour1
            colour= subperim(px(di(j)), py(di(j)));
            j= j+1;
        end;
        cutperim2x= px(di(j-1)) + kinkx(i) - bsize - 1;
        cutperim2y= py(di(j-1)) + kinky(i) - bsize - 1;
        
        % carry out cut 
        ri= regionprops(subcell,'solidity');
        if min([cutperim1x(1),cutperim1y(1),cutperim2x(1),cutperim2y(1)])>0 & ...
                max([cutperim1x(1),cutperim2x(1)])<=size(currcell,1) & ...
                max([cutperim1y(1),cutperim2y(1)])<=size(currcell,2)
            subcell= drawline(currcell, [cutperim1x(1) cutperim1y(1)],...
                [cutperim2x(1) cutperim2y(1)], 0);
            subcell= bwlabel(subcell, 4);
        else
            subcell=currcell;
        end
        
        % check cut
        rf= regionprops(subcell, 'majoraxislength', 'solidity');
        if min([rf.Solidity]) - min([ri.Solidity]) > 0.15
            % accept cut if it increases solidity
            currcell= subcell;   
        elseif min([rf.MajorAxisLength]) > mincelllength 
            % accept cut if it has not created too small cells
            currcell= subcell;
        else
            % ignore cut
            subcell= currcell;
        end;
        
    end;
    
end;

cutcell= zeros(size(cell));
cutcell(xmin:xmax, ymin:ymax)= subcell;    


