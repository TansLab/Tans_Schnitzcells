function [pos11]=DE_unfold_images(pos11,type)

%% DE: I didn't like that the grayscale figure was always beneath, and I
% wanted to put it side-by-side with the segmentation figure.
% INPUTS: 

% pos11 = [x1 y1 width height]
% position vector of newly created figure, here i interscept it and modify.

% type = 'gr'  for grayscale (will be on the right)
% type = 'seg' for the segmentation (will be on the left)

%%
% get figure's props:
width=pos11(3); height=pos11(4);
% get screen's props:
scren_size=get(0,'ScreenSize');
% calculate the offset for the left-bottom corner of the left figure:
sep=25;%separation
yoffset=(scren_size(4)-80-height)/2;
xoffset=(scren_size(3)-2*width-sep)/2;
% 80 pixels come into picture in y-coordinates because of the thick header
% of a window; sep = 25 (pixels) is the separation in x between the figures.

% coordinates of the bottom left corner:
bott_left_corner=[xoffset,scren_size(4)-yoffset-height-80]; 

switch type
    case 'seg'
        % new vector for pos11
        pos11=[bott_left_corner,width,height];
    case 'gr'
        bott_left_corner(1)=bott_left_corner(1)+width+sep;
        pos11=[bott_left_corner,width,height];
end