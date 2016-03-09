function cellcut = DJK_manual_cutOriginal(cell, cutx, cuty, figs);
% function cellcut = DJK_manual_cutOriginal(cell, cutx, cuty, oldcolor, newcolor);
%
% started from manual_kant, changed slightly
%
% breaks an individual cell by drawing line:
%

%--------------------------------------------------------------------------
% PARSE INPUT ARGUMENTS
%--------------------------------------------------------------------------
if nargin < 4
  figs = 0; % figs = 1 for graphical output (debug mode)
end

% if some error, will not cut, but will show figs
will_cut = 1; 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SETTINGS
%--------------------------------------------------------------------------
% box_radius = 4; % distance from clicked point where thin is used
% max_pixels_from_thin = 6.5; % will cut as far as so many pixels from thin
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------_
% GET SUBREGION CONTAINING CELL ONLY
%--------------------------------------------------------------------------
[fx,fy] = find(cell);
xmin = max(min(fx)-5,1);
xmax = min(max(fx)+5,size(cell,1));
ymin = max(min(fy)-5,1);
ymax = min(max(fy)+5,size(cell,2));
subcell = +cell(xmin:xmax, ymin:ymax);
subcutx = cutx-xmin+1;
subcuty = cuty-ymin+1;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% GET PERIMETER OF SUBCELL
%--------------------------------------------------------------------------
% perim is perimeter of dilated cell
subcell_dilated = imdilate(subcell,strel('disk',1));
perim = subcell_dilated & ~subcell;

% In Original: perim = bwperim(imdilate(subcell,strel('disk',1)));
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% FIND SMALLEST BOX AROUND CLICKED POINT WITH PERIMETERS OF BOTH SIDES
%--------------------------------------------------------------------------
% starting from clicked point, will increase a box untill 2 sides are
% found: this will be perims
perim_around_click = zeros(size(perim));
radius = 1;
while max2(perim_around_click)<2 & radius<30
  pxmin = max(cutx-xmin+1-radius,1);
  pxmax = min(cutx-xmin+1+radius,size(perim,1));
  pymin = max(cuty-ymin+1-radius,1);
  pymax = min(cuty-ymin+1+radius,size(perim,2));
  perim_around_click(pxmin:pxmax,pymin:pymax) = bwlabel(perim(pxmin:pxmax,pymin:pymax));
  radius = radius+1;
end

% % DJK: tried with circle, but slow
% % starting from clicked point, will increase a disk untill 2 sides are
% % found: this will be perims% mask_seed = zeros(size(subcell));
% mask_seed(round(subcutx), round(subcuty)) = 1;
% perim_around_click = zeros(size(subcell)); radius = 1;
% while max2(perim_around_click)<2 & radius<30
%   mask = imdilate(mask_seed, strel('disk',radius,0));
%   perim_around_click(find(mask>1)) = bwlabel(subcell(find(mask>1)));
%   radius = radius+1;
%   figure(21); imshowlabel(mask); pause; close(21);
% end

% if only 1 side is found, will not cut
if max2(perim_around_click)<2
  will_cut = 0; disp('perim_around_click with 2 sides not found');
end
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
% GET CLOSEST PIXEL ON PERIMETER AROUND CLICK
%-------------------------------------------------------------------------
[p1x,p1y] = find(perim_around_click>0);
if length(p1x)>1
  d1 = sqrt( (p1x - subcutx) .^ 2 + (p1y - subcuty) .^ 2);
  [ds1, dsI1]= sort(d1);
  px1 = p1x(dsI1(1));
  py1 = p1y(dsI1(1));
  color1 = perim_around_click(px1,py1);
else
  will_cut = 0;
end

%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% GET CLOSEST PIXEL ON PERIMETER AROUND CLICK ON OTHER SIDE
%-------------------------------------------------------------------------
if length(p1x)>1
  % remove pixels of side 1
  perim_around_click_1side = perim_around_click;
  perim_around_click_1side(perim_around_click_1side==color1) = 0;

  % get closest remaining pixel
  [p2x,p2y] = find(perim_around_click_1side>0);
  if length(p2x)>0
    d2 = sqrt( (p2x - subcutx) .^ 2 + (p2y - subcuty) .^ 2);
    [ds2, dsI2]= sort(d2);
    px2 = p2x(dsI2(1));
    py2 = p2y(dsI2(1));
  else
    will_cut = 0;
  end
end
%-------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PERFORM CUT
%--------------------------------------------------------------------------
if length(p1x)>1 & length(p2x)>1 
  % if click is more than 20 pixels to perimeter, don't cut
  if ds1(1) > 20 | ds2(1) > 20
    will_cut = 0;
  end
end

if (will_cut)
  % cut cell by drawing line
  subcellcut = drawline(subcell, [px1 py1],[px2 py2], 0); % disp('cut');
else
  subcellcut = subcell;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% ENLARGE OUTPUT
%--------------------------------------------------------------------------
% give unique colors
subcellcut = bwlabel(subcellcut, 4);

% enlarge
cellcut = zeros(size(cell));
cellcut(xmin:xmax,ymin:ymax) = subcellcut; 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% FIGURE OUTPUT
%--------------------------------------------------------------------------
if figs
  figure(111);

  % plot output cell
  subplot(2,1,2);
  PN_imshowlabel(p,subcellcut,[],[],[]);

  % plot original cell, perim, click point and line point
  subplot(2,1,1);
  output = subcell;
  output(find(perim>0)) = 2;
  output(find(perim_around_click>0)) = 3;
  output(round(subcutx),round(subcuty)) = 4;
  if length(p1x)>0
    output(round(px1),round(py1)) = 5; 
    if length(p2x)>0, output(round(px2),round(py2)) = 6; end
  end
  PN_imshowlabel(p,output,[],[],[]);

  pause; close(111);
end
%--------------------------------------------------------------------------
