% DJK_imagePlace writes image into other image
% only pixel values > 0 in source are written into target
% x,y indicates where center of im_source goes in im_target
%
% REQUIRED ARGUMENTS:
% 'im_target'   target image
% 'im_source'   source image
% 'x'           x location where source image will be put (note: x = col)
% 'y'           y location where source image will be put (note: y = row)

function im_output = DJK_imagePlace(im_target, im_source, x, y);

[rows1, cols1] = size(im_target);
[rows2, cols2] = size(im_source);    

row = round( y - 0.5*rows2 );
col = round( x - 0.5*cols2 );

% Find min row and column of im2 that will appear in im1
rmin2 = max(1,1-row);
cmin2 = max(1,1-col);    

% Find min row and column within im1 that im2 covers
rmin1 = max(1,1+row);
cmin1 = max(1,1+col);    

% Find max row and column of im2 that will appear in im1    
rmax2 = min(rows1-row, rows2);
cmax2 = min(cols1-col, cols2);    

% Find max row and column within im1 that im2 covers    
rmax1 = min(rows2+row, rows1);
cmax1 = min(cols2+col, cols1);    

% Simply copy im_target to im_output
im_output = im_target;  
 
% Check for the case where there is no overlap of the images
if rmax1 < 1 | cmax1 < 1 | rmax2 < 1 | cmax2 < 1 | rmin1 > rows1 | cmin1 > cols1 | rmin2 > rows2 | cmin2 > cols2
  return;
end

% Place im_number into im_target
im_temp = zeros(size(im_target));
im_temp(rmin1:rmax1, cmin1:cmax1) = im_source(rmin2:rmax2, cmin2:cmax2);
im_output(find(im_temp>0)) = im_temp(find(im_temp>0));