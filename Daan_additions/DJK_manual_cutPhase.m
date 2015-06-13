function cellcut = DJK_manual_cutPhase(cell, phcell, cutx, cuty, figs);
% function cellcut = DJK_manual_cutPhase(cell, phcell, cutx, cuty, oldcolor, newcolor);
%
% started from breakcell, but changed completely
%
% breaks an individual cell by approximately 2 pixels:
%  * first `thin' of the cell is taken close to click location
%  * phase image is checked along line between 2 endpoints (spurs) of this thin
%  * 2 pixels closest to maximum along this line are deleted 
%  * this is repeated for parallel lines untill cell is cut
%

%--------------------------------------------------------------------------
% PARSE INPUT ARGUMENTS
%--------------------------------------------------------------------------
if nargin == 4
  figs = 0; % figs = 1 for graphical output (debug mode)
end

% make sure working with double cell containing 1s and 0s
cell = +(cell > 0);

% check euler number
r = regionprops(cell, 'eulernumber');
if [r.EulerNumber] ~= 1
  disp(['Circular region in thin, might give problems. Euler Number= ',num2str([r.EulerNumber])]);
end

% if some error, will not cut, but will show figs
will_cut = 1; 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SETTINGS
%--------------------------------------------------------------------------
box_radius = 4; % distance from clicked point where thin is used
max_pixels_from_thin = 6.5; % will cut as far as so many pixels from thin
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% GET SUBREGION CONTAINING CELL ONLY
%--------------------------------------------------------------------------
[fx,fy] = find(cell);
xmin = max(min(fx)-5,1);
xmax = min(max(fx)+5,size(cell,1));
ymin = max(min(fy)-5,1);
ymax = min(max(fy)+5,size(cell,2));
subcell = cell(xmin:xmax, ymin:ymax);
subphcell = phcell(xmin:xmax, ymin:ymax);
subcutx = cutx-xmin+1;
subcuty = cuty-ymin+1;
%figure(111);imshowlabel(subcell);pause;close(111);
%figure(111);imshow(imadjust(subphcell));pause;close(111);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% FIND THIN (CENTRE LINE) OF CELL, ALSO CLOSE TO CLICK POINT
%--------------------------------------------------------------------------
subthin = bwmorphmelow(subcell, 'thin', inf); % make thin of image
%figure(111);imshowlabel(subthin);pause;close(111);
subthin = bwmorphmelow(bwmorphmelow(subthin, 'diag', 1), 'thin', inf); % clean up thin (remove spurious spurs)
%figure(111);imshowlabel(subthin);pause;close(111);

% ONLY TAKE REGION AROUND CLICKED POINT
% 2DO make it a circle instead of cube
xmin_box = round( max( subcutx - box_radius, 1              ) );
xmax_box = round( min( subcutx + box_radius, size(subcell,1)) );
ymin_box = round( max( subcuty - box_radius, 1              ) );
ymax_box = round( min( subcuty + box_radius, size(subcell,2)) );
subthin_box = zeros(size(subthin));
subthin_box(xmin_box:xmax_box,ymin_box:ymax_box) = subthin(xmin_box:xmax_box,ymin_box:ymax_box);
%figure(111);imshowlabel(subthin_box);pause;close(111);

% DETERMINE ANGLE OF subthin_box
r = regionprops(subthin_box, 'Orientation');
if length(r)~=0 % to avoid error when r is empty, when clicked at end of cell
  angle = r.Orientation;
  angle_perp = angle - 90; % 90 - abs(angle);
  angle_radians = ((2*pi)/360)*angle;
  angle_perp_radians = ((2*pi)/360)*angle_perp;
  offset_step_x = -sin(angle_perp_radians);
  offset_step_y =  cos(angle_perp_radians);
  if (figs)
    %disp(['angle and angle_perp are : ' num2str(angle) ' ' num2str(angle_perp)]);
    %disp(['offset_step_x and y are : ' num2str(offset_step_x) ' ' num2str(offset_step_y)]);
  end
end

% FIND SPUR POINTS (END POINTS OF SUBTHIN_BOX, WILL DEFINE LINE)
spurs = subthin_box & ~bwmorphmelow(subthin_box, 'spur', 1);
[sx, sy] = find(spurs > 0);
if length(sx) ~= 2 
  disp(['More or less than 2 spur points = ',num2str(length(sx))]);
  disp(['Not cutting anything. Try again or change box_radius.']);
  will_cut = 0;
end
%--------------------------------------------------------------------------


% 2DO, zelf aanklikken tussen welke punten gefit moet worden
% http://blogs.mathworks.com/pick/2007/12/26/advanced-matlab-buttondownfcn/


%--------------------------------------------------------------------------
% CHECK PHASE CONTRAST VALUES ALONG LINE BETWEEN 2 SPUR POINTS AND
% DETERMINE CUT POINT
%--------------------------------------------------------------------------
if will_cut
  % get phase data along line
  avp = improfile(subphcell, sy, sx , 4*box_radius, 'bicubic'); %bilinear is faster
  avp = avp'/median(avp); %normalize

% if isnan(avp) % probably clicked at end of cell, so no good thin
%   disp(['Could not get improfile, probably no good thin. Not cutting.']);
%   will_cut = 0;
%   fitCoef = 0; max_fit = 0; max_ind = 0;
% else
  % fit parabola through center of data and determine maximum
  fitx = [(box_radius+1):(3*box_radius)]; % center half of data
  fitCoef = polyfit(fitx, avp(fitx), 2); % fit parabola
  max_fit = -fitCoef(2)/(2*fitCoef(1)); % maximum of parabola

  % if fit is bad, do not cut
  if (fitCoef(1)>0)
    disp(['Not cutting, cause parabola fit to ph contrast is not nice']);
    will_cut = 0;
  end

  % if max_fit is far away from clicked point, do not cut
  if (max_fit<box_radius+1 | max_fit>3*box_radius)
    disp(['Not cutting, cause max_fit is far away from clicked point']);
    will_cut = 0;
  end

  % find local maxima in phase (potential point to be cut)
  [max_val, max_ind] = max(avp);

  % determine coordinates of maximum along this line 
  % 2DO use wider line
  cut_point_x  = sx(1) + ((max_fit -1)/(4*box_radius))*(sx(2)-sx(1) );
  cut_point_y  = sy(1) + ((max_fit -1)/(4*box_radius))*(sy(2)-sy(1) );
  cut_point_xx(1) = cut_point_x - 0.5*sin(angle_radians);
  cut_point_yy(1) = cut_point_y + 0.5*cos(angle_radians);
  cut_point_xx(2) = cut_point_x + 0.5*sin(angle_radians);
  cut_point_yy(2) = cut_point_y - 0.5*cos(angle_radians);
%   cut_point_x3 = sx(1) + ((max_fit-1)/(4*box_radius))*(sx(2)-sx(1) ); % old way
%   cut_point_y3 = sy(1) + ((max_fit-1)/(4*box_radius))*(sy(2)-sy(1) ); % old way
%   cut_point_x4 = sx(1) + ((max_fit+1)/(4*box_radius))*(sx(2)-sx(1) ); % old way
%   cut_point_y4 = sy(1) + ((max_fit+1)/(4*box_radius))*(sy(2)-sy(1) ); % old way
  % figure(111);hold on;
  % plot(cut_point_y, cut_point_x, 'ro');
  % plot([cut_point_yy(1) cut_point_yy(2)], [cut_point_xx(1) cut_point_xx(2)], 'bo');
  % plot([cut_point_y3 cut_point_y4], [cut_point_x3 cut_point_x4], 'go');
  % set(gca,'YDir','reverse')
  % pause;close(111);


  % % perpendicular line through this point will have angle 90 degress - angle
  % angle_perp = 90 - angle;
else
  avp = []; fitCoef = 0; max_fit = 0; % otherwise fig error
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PERFORM CUT
%--------------------------------------------------------------------------
% used to check whether cutpts are in cell
subcell_idx = find(subcell);

% subperim is perimeter of cell (will cut untill perimeter is reached)
subperim = imdilate(subcell,strel('disk',1)) & ~subcell; %perim = bwperim(subcell); % perim = bwperim(imdilate(subcell, strel('disk',1)));
subperim_idx = find(subperim);
%figure(112);imshowlabel(perim);pause;

% cut_pts will contain pixels that will be cut
cut_pts = zeros(size(subthin));

if (will_cut)
  for i=[1:4] % 1: xx(2)-offset  2: xx(1)-offset 3: xx(2)+offset  4: xx(1)+offset
    %disp(['i ' num2str(i)]);
    perim_reached = false; offset = 0;
    while (~perim_reached & abs(offset) < max_pixels_from_thin) % loop untill perim is reached

      % offset pixels
      offset_x = cut_point_xx(1+rem(i,2)) + offset*offset_step_x;
      offset_y = cut_point_yy(1+rem(i,2)) + offset*offset_step_y;
      %disp(['(' num2str(offset_x) ',' num2str(offset_y) ')']);

      % round of pixels
      offset_x = round(offset_x); offset_y = round(offset_y);
      %disp(['(' num2str(offset_x) ',' num2str(offset_y) ')']);

      if ( offset_x>0 & offset_y>0 & offset_x<=size(cut_pts,1) & offset_y<=size(cut_pts,2) )
        % put pixels in cut_pts if they correspond to a location in subcell
        if ( find(subcell_idx==sub2ind(size(subcell), offset_x, offset_y)) > 0)
          cut_pts(offset_x,offset_y) = 1;
        end
        %figure(111);imshowlabel(cut_pts);pause;close(111);

        % check whether perim is reached
        if (find(subperim_idx==sub2ind(size(subcell), offset_x, offset_y)) > 0)
          perim_reached = true;
        end
      else
        perim_reached = true;
      end

      if (i>2)
        offset = offset + 0.5;
      else 
        offset = offset - 0.5;
      end

      %figure(111);imshowlabel(bwlabel(subcell & ~cut_pts, 4));pause;close(111);
    end
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% ENLARGE OUTPUT
%--------------------------------------------------------------------------
% remove cut points
subcellcut = subcell;
subcellcut(find(cut_pts>0)) = 0;

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
  figure(111); clf;
  subplot(2,1,1);
  x = 1:length(avp);
  plot(x, avp, 'bo', x, polyval(fitCoef,x),'k-', [max_fit, max_fit], [0.9, 1.1], 'r-'); %, [max_ind, max_ind], [0.9, 1.1], 'g-'
  legend('phase','poly fit','max fit',-1);
  xlim([0 length(avp)+1]);
  ylim([0.6 1.4]);

  % plot phase with thin and thin2
  subplot(2,2,3); 
  [subphcell_index, map_gray] = gray2ind(imadjust(subphcell), 64);
  map_size = size(map_gray);
  map_gray(map_size(1)+1,:) = [ 0 0 1 ];
  map_gray(map_size(1)+2,:) = [ 1 0 0 ];
  map_gray(map_size(1)+3,:) = [ 0 1 0  ];
  subphcell_index(find(subthin>0)) = map_size(1);
  subphcell_index(find(cut_pts>0)) = map_size(1)+2;
  subphcell_index(find(subthin_box>0)) = map_size(1)+1;
  subimage(subphcell_index, map_gray);
  axis off;

  % plot output cell
  subplot(2,2,4);
  PN_imshowlabel(p,subcellcut,[],[],[]);

  pause;close(111);
end
%--------------------------------------------------------------------------
