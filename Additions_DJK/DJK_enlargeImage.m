function im_out = DJK_enlargeImage(im_in, fullSize, rect, stabilize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enlargeImage puts image back in original size
% if stabilize, will center in original size (less shaking)
% rect is [y_top, x_left, y_bottom, x_right] , x and y hier door elkaar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
or_size = size(im_in);
x = round( rect(2) + or_size(2)/2 );
y = round( rect(1) + or_size(1)/2 );
im_out = DJK_imagePlace(zeros(fullSize), im_in, x, y);

if (stabilize) 
  [fy, fx]= find(im_out>0);
  x = round( fullSize(2) - mean(fx) );
  y = round( fullSize(1) - mean(fy) );
  im_out = DJK_imagePlace(zeros(fullSize), im_out, x, y);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before 090414 used stabilizing based purely on rect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% im_out = zeros(fullSize);
% or_size = size(im_in);
% off_set_x = rect(1);
% off_set_y = rect(2);
% if (stabilize) 
%   off_set_x = uint16( (0.5+fullSize(1)+rect(1)-rect(3))/2 );
%   off_set_y = uint16( (0.5+fullSize(2)+rect(2)-rect(4))/2 );
% end
% for x = 1:or_size(1)
%     for y = 1:or_size(2)
%         im_out(x+off_set_x, y+off_set_y) = im_in(x,y);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






