
function cutcell= breakcell(cell, phcell, maxthresh, minthresh, mincelllength, figs, thinx, thiny)
% function cutcell= breakcell(cell, phcell, maxthresh, minthresh, mincelllength, figs, thinx, thiny)
%
% breaks an individual cell by running a box along the `thin' of cell
% cell is cut at points where there is a minimum in the segmented image (cell) 
%  or a maximum in the phase image (phcell)
% maxthresh & minthresh:    threshold values (smaller they are the more cuts are included)
% mincelllength:    cuts which create cells smaller than this will be ignored
% 
% rather than calculate thin from cell, can receive coordinates of thin via (thinx, thiny) 

if nargin == 5
    figs= 0; % figs= 1 for graphical output (debug mode)
end;
maxdist= 4; % distance from pt used to test if it is a maximum
mindist= 4; % distance from pt used to test if it is a minimum
diffthresh= 3; % size of jump in diff of minpts or maxpts 

% make sure image is black and white and not logical
cellno= max2(cell);
cell= +(cell > 0);
phcell= +phcell;

% check Euler Number
if nargin < 7
    r= regionprops(cell, 'eulernumber');
    if [r.EulerNumber] ~= 1
        if figs==1
            disp(['Circular region in thin. Euler Number= ',num2str([r.EulerNumber])]);
        end;
        cutcell= bwlabel(cell, 4);
        return;
    end;
end;

% extract subimages
[fx, fy]= find(cell);
extra= 5;
xmin= max(min(fx) - extra, 1);
xmax= min(max(fx) + extra, size(cell,1));
ymin= max(min(fy) - extra, 1);
ymax= min(max(fy) + extra, size(cell,2));
subcell= cell(xmin:xmax, ymin:ymax);
originalsubcell= subcell;
subphcell= phcell(xmin:xmax, ymin:ymax);


% FIND THIN (CENTRE LINE) OF CELL

if nargin < 7
    
    % make thin of image
    thin= bwmorphmelow(subcell, 'thin', inf);
    % clean up thin (remove spurious spurs)
    thin= bwmorphmelow(bwmorphmelow(thin, 'diag', 1), 'thin', inf);
    
    % find spur points
    spurs= thin & ~bwmorphmelow(thin, 'spur', 1);
    [sx, sy]= find(spurs > 0);
    if figs == 1
        disp(['No of spur points in thin= ',num2str(length(sx))]);
    end;
    if length(sx) > 2
        for i= 2:length(sx)
            subcell2= zeros(size(subcell));
            subcell2(sx(1), sy(1))= 1;
            subcell2(sx(i), sy(i))= 1;
            subcell2= imdilate(subcell2, strel('disk', 2));
            subcell3= subcell2 | thin;
            subcell4= bwmorphmelow(subcell3, 'spur', inf);
            subcell5= bwmorphmelow(subcell4, 'thin', inf);
            subcellthin{i-1}= bwmorphmelow(subcell5, 'spur', 2);
        end;
    else
        subcellthin{1}= thin;
    end;
    
else
    
    % create thin image
    thin= zeros(size(subcell));
    fx= thinx - xmin + 1;
    fy= thiny - ymin + 1;
    for i= 1:length(fx)
        thin(fx(i), fy(i))= 1;
    end;
    subcellthin{1}= thin;
    
    % find spur points
    spurs= thin & ~bwmorphmelow(thin, 'spur', 1);
    [sx, sy]= find(spurs > 0);
    if length(sx) ~= 2
        disp('Given thin has multiple spur points.');
        cutcell= cell;
        return;
    end;
    
end;

% RUN THROUGH EACH THIN TO FIND POINTS TO CUT

for k= 1:length(sx)-1
    
    % obtain points ordered along thin of image
    clear avp avs;
    if nargin < 7
        [fx, fy]= walkthin(subcellthin{k});
        if isempty(fx)
            continue;
        end;
    end;
    j= 1;
    for bsize= 2:5
        for i= 1:length(fx),
            % extract mean of pixels in phase image lying in box bsize
            axmin= max(1, fx(i) - bsize);
            axmax= min(size(subcell,1), fx(i) + bsize);
            aymin= max(1, fy(i) - bsize);
            aymax= min(size(subcell,2), fy(i) + bsize);
            avp(j,i)= mean2(subphcell(axmin:axmax, aymin:aymax));
            avs(j,i)= mean2(originalsubcell(axmin:axmax, aymin:aymax));
        end;
        j= j+1;
    end;
    % normalize 
    avs= mean(avs)/median(median(avs));
    avp= mean(avp)/median(median(avp));
    
    
    % find local maxima in phase (potential points to be cut)
    maxpts= [];
    for i= 1 + maxdist:length(avp) - maxdist
        if (avp(i - maxdist) < avp(i)) & (avp(i + maxdist) < avp(i))
            maxpts= [maxpts i];
        end;
    end; 
    if ~isempty(maxpts)
        % maxpts contains multiple points for each maximum
        % find boundaries between sets of points
        maxboundaries= unique([1 find(diff(maxpts) > diffthresh) length(maxpts)]);
        
        if length(maxpts) == 1
            maxs{k}= maxpts;
            % measure steepness of maximum
            maxsc{k}= mean([avp(maxs{k}) - avp(maxs{k}-maxdist),...
                    avp(maxs{k}) - avp(maxs{k}+maxdist)]);
        else
            % each maximum is the average of each set of points associated with it
            for i= 2:length(maxboundaries)
                submaxpts= maxpts(maxboundaries(i-1) + 1:maxboundaries(i));
                [av2m, av2mi]= max(avp(submaxpts));
                maxs{k}(i-1)= submaxpts(av2mi);
                % measure steepness of maximum
                maxsc{k}(i-1)= mean([avp(maxs{k}(i-1)) - avp(maxs{k}(i-1)-maxdist),...
                        avp(maxs{k}(i-1)) - avp(maxs{k}(i-1)+maxdist)]);
            end;
        end;
        
        % choose steep maxima only
        cutptsmax{k}= [maxs{k}(find(maxsc{k} > maxthresh))];
    else
        cutptsmax{k}= [];
        maxsc{k}= [];
    end;
    
    
    
    
    % find local minima in segmented (potential points to be cut)
    minpts= [];
    for i= 1 + mindist:length(avs) - mindist
        if (avs(i - mindist) > avs(i)) & (avs(i + mindist) > avs(i))
            minpts= [minpts i];
        end;
    end; 
    if ~isempty(minpts)
        % minpts contains multiple points for each minimum
        % find boundaries between sets of points
        minboundaries= unique([1 find(diff(minpts) > diffthresh) length(minpts)]);
        
        if length(minpts) == 1
            mins{k}= minpts;
            % measure steepness of minimum
            minsc{k}= mean([avs(mins{k}-mindist) - avs(mins{k}),...
                    avs(mins{k}+mindist) - avs(mins{k})]);
        else
            % each minimum is the average of each set of points associated with it
            for i= 2:length(minboundaries)
                subminpts= minpts(minboundaries(i-1) + 1:minboundaries(i));
                [av2m, av2mi]= min(avs(subminpts));
                mins{k}(i-1)= subminpts(av2mi);
                % measure steepness of minimum
                minsc{k}(i-1)= mean([avs(mins{k}(i-1)-mindist) - avs(mins{k}(i-1)),...
                        avs(mins{k}(i-1)+mindist) - avs(mins{k}(i-1))]);
            end;
        end;
        
        
        % choose large minima only
        cutptsmin{k}= [mins{k}(find(minsc{k} > minthresh))];
    else
        cutptsmin{k}= [];    
        minsc{k}= [];
    end;
   
    
   % find cutpts
   cutpts{k}= unique([cutptsmin{k} cutptsmax{k}]);

    
    
    % CUT CELLS
    
    if ~isempty(cutpts{k})
        % must cut cells
        
        % sort cutpts in order of distance from centre of thin
        [cs, csi]= sort(abs(cutpts{k} - round(length(fx)/2)));
        cutpts{k}= cutpts{k}(csi);
        
        cutx= fx(cutpts{k});
        cuty= fy(cutpts{k});
        
        % now divide the cell into seperate cells by cutting it across
        perim= bwperim(imdilate(subcell, strel('disk',1)));
       
        for i= 1:length(cutpts{k})
            
            bsize= 8;
            sxmin= max(1, cutx(i) - bsize);
            sxmax= min(size(perim,1), cutx(i) + bsize);
            symin= max(1, cuty(i) - bsize);
            symax= min(size(perim,2), cuty(i) + bsize);
            subperim= perim(sxmin:sxmax, symin:symax);
            [subperim, noperims]= bwlabel(subperim);
            
            % if noperims ~= 1  % JCR: noperims==0 causes problems below
            if noperims > 1
                
                % cutpt is not near end of cell. Go ahead and cut.
                currcell= subcell;
                [px, py]= find(subperim> 0);
                
                % find distances to perimeter from cutpt
                d= sqrt((px - bsize - 1).^2 + (py - bsize - 1).^2);
                [ds, di]= sort(d);
                
                % find first cutting point on perimeter
                cutperim1x= px(di(1)) + sxmin - 1;
                cutperim1y= py(di(1)) + symin - 1;
                colour1= subperim(px(di(1)), py(di(1)));
                
                % find second cutting point on perimeter
                colour= colour1;
                j= 2;
                while colour == colour1
                    colour= subperim(px(di(j)), py(di(j)));
                    j= j+1;
                end;
                cutperim2x= px(di(j-1)) + sxmin - 1;
                cutperim2y= py(di(j-1)) + symin - 1;
                
                % carry out cut 
                subcell= drawline(currcell, [cutperim1x(1) cutperim1y(1)],...
                    [cutperim2x(1) cutperim2y(1)], 0);
                subcell= bwlabel(subcell, 4);
                
                % check cut
                rf= regionprops(subcell, 'majoraxislength');
                if min([rf.MajorAxisLength]) > mincelllength 
                    % accept cut if it has not created too small cells
                    currcell= subcell;
                else
                    % ignore cut
                    subcell= currcell;
                end;
                
            end;
        end;   
    else
        
        cutx= [];
        
    end;
    
end;
            
cutcell= zeros(size(cell));
cutcell(xmin:xmax, ymin:ymax)= subcell;    
    
    
% OUTPUT FIGURES

% figure output
if figs == 1
    
    figure; clf; 
    title(cellno);
    
    subplot(3,1,1);
    plot(1:length(avp), avp, 'b.', 1:length(avs), avs, 'k.');
    legend('phase','segmented',-1);
    xlabel(['cell number ', num2str(cellno)]);
    
    % plot phase and thin
    subplot(3,1,2);
    cutptsim= zeros(size(thin));
    if ~isempty(cutx)
        for i= 1:length(cutx)
            cutptsim(cutx(i), cuty(i))= 1;
        end
        imshow(makergb2(subphcell, thin, cutptsim));
    else 
        imshow(makergb2(subphcell, thin));
    end;

    % plot output cell
    subplot(3,1,3);
    imshowlabel(subcell);
    
    for i= 1:length(sx)-1
        disp([num2str(i),': maxsc(',num2str(size(maxsc{i},2)),') ',num2str(maxsc{i})]);
        disp([num2str(i),': minsc(',num2str(size(minsc{i},2)),') ',num2str(minsc{i})]);
    end;
end;