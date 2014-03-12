function cutcell= DJK_cutcurvbpts(cell, maxthresh, maxcellwidth, figs)
% function cutcell= DJK_cutcurvbpts(cell, maxthresh, maxcellwidth, figs)
%
% copied from cutcurvbpts.
%
% cuts cells at points where both sides of the cell are sufficiently concave
%   maxthresh:      smaller maxtresh the more points are cut
%   maxcellwidth:   two points must be closer than this for cut to be accepted
%
% finds branchpoints in cell (points with more than 2 neighbours on thin)
%   and extracts and cuts subregions, one for each branchpoint
% hence can handle multiply connected cells

if nargin == 3
    figs= 0;
end;
maxdist= 4; % distance from pt used to test if it is a maximum
diffthresh= 3; % size of jump in diff of minpts or maxpts 
holesize= 20;

% DJK: 090111 uitgezet, want teveel fill in 081011 pos3crop 210,212:218
% % fill in any holes in cell
% cc= bwlabel(~cell, 4);
% r= regionprops(cc, 'area');
% sinpixel= find([r.Area] < holesize);
% for i= 1:length(sinpixel)
%     [fx, fy]= find(cc == sinpixel(i));
%     for j= 1:length(fx)
%         cell(fx(j), fy(j))= 1;
%     end;
% end;

% DJK 090310 uitgezet, want gets stuck in 2009-03-03 pos4 improved_edge123_fill3_pix275_hat6 seg008
%cell= +(rmsinglepointconnections(cell) > 0);
currcell= cell;

% FIND BRANCHPOINTS
thin= bwmorphmelow(cell, 'thin', inf);
filt = ones(3,3);
filt(2,2) = 0;

bpts= imfilter(thin, filt);
bpts= (bpts > 2) & thin;
bpts= bwmorphmelow(bpts, 'shrink', inf);

[bx, by]= find(bpts > 0);

% RUN THROUGH BRANCHPOINTS
for nb= 1:length(bx)
    %disp(['   branch point ',num2str(nb)]);
    
    extra= 15;
    ex= 0;
    % check that perim does not have any ends
    % (these could arise from badly cutoff cells -> change extra)
    while ~isempty(ex) & extra < 20
        % EXTRACT SUBIMAGE
        clear subcell;
        
        
        xmin = min(bx(nb))-extra;
        xmax= min(max(bx(nb)) + extra, size(cell,1));
        ymin = min(by(nb))-extra;
        ymax= min(max(by(nb)) + extra, size(cell,2));
        
        locx = extra+1;
        locy = extra+1;
        
        if xmin < 1,
%           locx = locx + (1-xmin);  JCR: is this a bug?
            locx = locx - (1-xmin);
            xmin = 1;
        end;
        if ymin < 1,
%           locy = locy + (1-ymin);  JCR: is this a bug?
            locy = locy - (1-ymin);
            ymin = 1;
        end;
        
        %         xmin= max(min(bx(nb)) - extra, 1);
        %         ymin= max(min(by(nb)) - extra, 1);
        subcell= bwlabel(cell(xmin:xmax, ymin:ymax), 4);
        %       [delx, dely]= find(subcell ~= 0 & subcell ~= subcell(extra, extra)); 
        % nitzan addition: change extra to extra+1 to check point bx,by and not to miss it by one pixel.
      
        % added by ME on 2002-07-27, because sometimes the subcell i smaller
        % than locx or locy...
        if locx <= size(subcell,1) & locy <= size(subcell,2),
            [delx, dely]= find(subcell ~= 0 & subcell ~= subcell(locx, locy));
        else
            cutcell = currcell;
            disp('WARNING: early return from cutcurvbpts!');
            return
        end;
       
        for i= 1:length(delx)
            subcell(delx(i), dely(i))= 0;
        end;
        subcell= +(subcell > 0);
        
        %added by ME on 2002-12-03 due to weird problems:
              if max2(subcell)==0,
                      cutcell = currcell;
                      disp('WARNING: early return from cutcurvbpts!');
                      return;
              end;
        
        % FIND PERIMETER
        [px, py]= tracecontour(subcell);
        perimp= zeros(size(subcell));
        for i= 1:length(px)
            perimp(px(i), py(i))= 1;
        end;
        endpts= imfilter(perimp, filt);
        [ex, ey]= find(endpts == 1 & perimp);
        extra= extra + 1;
    end;

    % FIND INTERNAL ANGLES
    inc= 5;
    clear a;
    a=0;
% nitzan addition - sometimes px is too short, then a is undefined and routine gets stuck in line 135.
    if length(px)<3*inc inc=(length(px)-mod(length(px),3))/3;end
% end nitzan addition.
    for i= inc+1:length(px)-inc
        
        % find forward and backward vectors
        v1= [px(i+inc) - px(i) ; py(i+inc) - py(i)];
        v2= [px(i-inc) - px(i) ; py(i-inc) - py(i)];
        if sqrt(sum(v1.^2))>0 v1= v1/sqrt(sum(v1.^2));end
        if sqrt(sum(v2.^2))>0 v2= v2/sqrt(sum(v2.^2));end
        
        % find angle between them
        a(i)= real(acos(sum(v1.*v2)));
        
        % find centre vector
        vm= [(px(i+inc) + px(i-inc))/2 - px(i) ; (py(i+inc) + py(i-inc))/2 - py(i)];
        lenvm= sqrt(sum(vm.^2));
        if lenvm > 0
            
            vm= round(2*vm/lenvm);
            ptx= [px(i) - vm(1) ; px(i) + vm(1)];
            pty= [py(i) - vm(2) ; py(i) + vm(2)];
            sc= sqrt((ptx - px(i+inc)).^2 + (pty - py(i+inc)).^2 ...
                + (ptx - px(i-inc)).^2 + (pty - py(i-inc)).^2);
            
            % determine if centre vector lies inside or outside of cell
            cx= ptx(minI(sc));
            cy= pty(minI(sc));
            if cx < 1 | cx > size(subcell, 1)
                sgn= 1;
            elseif cy < 1 | cy > size(subcell, 2)
                sgn= 1;
            else
                sgn= subcell(ptx(minI(sc)), pty(minI(sc)));
            end;
            % assign sign of angle
            if sgn == 0
                a(i)= -a(i);
            end;
            
        end;
    end;
    
    % POTENTIAL CUTS ARE FOR NEGATIVE ANGLES
    ineg= find(a < 0);
    b= a(ineg);
    
    % find local maxima (potential points to be cut)
    maxpts= [];
    for i= 1 + maxdist:length(b) - maxdist
        if (b(i - maxdist) < b(i)) & (b(i + maxdist) < b(i))
            maxpts= [maxpts i];
        end;
    end; 
    if ~isempty(maxpts)
        % maxpts contains multiple points for each maximum
        % find boundaries between sets of points
        maxboundaries= unique([1 find(diff(maxpts) > diffthresh) length(maxpts)]);
        
        if length(maxpts) == 1
            maxs= maxpts;
            % measure steepness of maximum
            maxsc= mean([b(maxs) - b(maxs-maxdist),...
                    b(maxs) - b(maxs+maxdist)]);
        else
            % each maximum is the average of each set of points associated with it
            maxs= []; maxsc= [];
            for i= 2:length(maxboundaries)
                submaxpts= maxpts(maxboundaries(i-1) + 1:maxboundaries(i));
                maxs(i-1)= submaxpts(maxI(b(submaxpts)));
                % measure steepness of maximum
                maxsc(i-1)= mean([b(maxs(i-1)) - b(maxs(i-1)-maxdist),...
                        b(maxs(i-1)) - b(maxs(i-1)+maxdist)]);
            end;
        end;
        
        % choose steep maxima only
        bcutpts= [maxs(find(maxsc > maxthresh))];
        cutpts= ineg(bcutpts); 
    else
        cutpts= [];
        maxsc= [];
    end;
    
    % PERFORM CUT
    tempcell= subcell;
    if length(cutpts) > 1
        
        % find cut point opposite each cut point
        for i= 1:length(cutpts)
            d= sqrt((px(cutpts(i)) - px(cutpts)).^2 + (py(cutpts(i)) - py(cutpts)).^2);
            [ds, dsI]= sort(d);
            % cutpts only accepted if they're close enough
            if ds(2) < maxcellwidth
                cutvec(i)= cutpts(dsI(2));
            else 
                cutvec(i)= 0;
            end;
        end;
        
        for i= 1:length(cutpts)
            if cutvec(i) ~= 0
                % check cutting points are reciprocal
                if cutpts(i) == cutvec(find(cutpts == cutvec(i)))
                    % do the cut
                    pt1 = [px(cutpts(i)) + xmin - 1, py(cutpts(i)) + ymin - 1];
                    pt2 = [px(cutvec(i)) + xmin - 1, py(cutvec(i)) + ymin - 1];
                    currcell= drawline(currcell, pt1, pt2, 0);
                    if figs == 1
                        tempcell= drawline(tempcell, [px(cutpts(i)), py(cutpts(i))], [px(cutvec(i)), py(cutvec(i))], 0);
                        pt_middle = (pt1+pt2)/2;
                        if nb == 10
                          disp(['pt1, pt2, pt_middle: (' num2str(pt1(1)) ',' num2str(pt1(2)) ') (' num2str(pt2(1)) ',' num2str(pt2(2)) ') (' num2str(pt_middle(1)) ',' num2str(pt_middle(2)) ')']);
                        end
                    end;
                end;
            end;
        end;
    end;
    
    if figs == 1 & nb == 10
        %clf;
        figure(50+nb);
        subplot(1,2,1);
        imshowlabel(bwlabel(tempcell, 4));

  
%         imshow(subcell);
%         hold on;plot(py, px, 'r.'); cutpts_x = px(cutpts); cutpts_y = py(cutpts); plot(cutpts_y, cutpts_x, 'w.');hold off;
  
        subplot(1,2,2);
        for i=1:length(px), subcell(px(i),py(i)) = 2; end
        cutpts_x = px(cutpts); cutpts_y = py(cutpts); 
        for i=1:length(cutpts_x), subcell(cutpts_x(i),cutpts_y(i)) = 3; end
        imshowlabel(subcell);
        %pause(5);
  
    end;

end;

cutcell= bwlabel(currcell, 4);