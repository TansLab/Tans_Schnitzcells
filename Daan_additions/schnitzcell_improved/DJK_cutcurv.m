function cutcell= DJK_cutcurv(cell, maxthresh, maxcellwidth, mincelllength, figs)
% function cutcell= DJK_cutcurv(cell, maxthresh, maxcellwidth, mincelllength, figs)
%
% copied from cutcurv, only changed figs part
%
% cuts cells at points where both sides of the cell are sufficiently concave
%   maxthresh:      smaller maxtresh the more points are cut
%   maxcellwidth:   two points must be closer than this for cut to be accepted
%   mincelllength:  cut ignored if it creates `cell' smaller than this

if nargin == 4
    figs= 0;
end;
maxdist= 4; % distance from pt used to test if it is a maximum
diffthresh= 3; % size of jump in diff of minpts or maxpts 
holesize= 20;
minperimeterlength= 12;

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

% cell= rmsinglepointconnections(cell);
currcell= cell;

% FIND PERIMETER
% perim= bwperim(cell);
% perim2= perim;
[px, py]= tracecontour(cell);
if length(px) < minperimeterlength
    cutcell= cell;
    return
end;


% FIND INTERNAL ANGLES
inc= 5;
clear a;
% nitzan addition, to take care of small cells. also: "end" at end of file.
if length(px)<3*inc
    cutcell=zeros(size(cell));
else
% end of addition, except the "end" at the end of the file.    
for i= inc+1:length(px)-inc
    
    % find forward and backward vectors
    v1= [px(i+inc) - px(i) ; py(i+inc) - py(i)];
    v2= [px(i-inc) - px(i) ; py(i-inc) - py(i)];
    v1= v1/sqrt(sum(v1.^2));
    v2= v2/sqrt(sum(v2.^2));
    
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
        if cx < 1 | cx > size(cell, 1)
            sgn= 1;
        elseif cy < 1 | cy > size(cell, 2)
            sgn= 1;
        else
            sgn= cell(ptx(minI(sc)), pty(minI(sc)));
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
                if ((px(cutpts(i))-px(cutvec(i)))^2+(py(cutpts(i))-py(cutvec(i)))^2)
                    % do the cut
                    currcell2= drawline(currcell, [px(cutpts(i)), py(cutpts(i))], ...
                        [px(cutvec(i)), py(cutvec(i))], 0);          
                    % nitzan addition - I moved this "end" to line 167.
                    %           end;
                    
                    % check cut
                    rf= regionprops(bwlabel(currcell2, 4), 'majoraxislength');
                    if min([rf.MajorAxisLength]) > mincelllength 
                        % accept cut if it has not created too small cells
                        currcell= currcell2;
                    end;
                    % nitzan addition - I moved the "end" to here from line 158.     
                end
            end;
        end;
    end;
end;

if figs == 1
  figure(51);
  subplot(1,2,1);
  imshowlabel(bwlabel(currcell, 4));

  subplot(1,2,2);  %figure(52); 
  cell = +cell;
  for i=1:length(px), cell(px(i),py(i)) = 2; end
  cutpts_x = px(cutpts); cutpts_y = py(cutpts); 
  for i=1:length(cutpts_x), cell(cutpts_x(i),cutpts_y(i)) = 3; end
  imshowlabel(cell); 

  for i=1:length(cutpts_x)
    disp(['cutpts: (' num2str(cutpts_x(i)) ',' num2str(cutpts_y(i)) ')']);
  end
  
  %keyboard; %pause; close(51); %close(52);
%   [bwperim_x, bwperim_y] = find(bwperim(cell)>0);
%   for i=[1:length(px)]
%     if i>length(bwperim_x)
%       disp([ '(' str3(0) ',' str3(0) ') vs (' str3(px(i)) ',' str3(py(i)) ')']); 
%     else
%       disp([ '(' str3(bwperim_x(i)) ',' str3(bwperim_y(i)) ') vs (' str3(px(i)) ',' str3(py(i)) ')']); 
%     end
%   end
end;


cutcell= bwlabel(currcell, 4);

% nitzan addition:
end