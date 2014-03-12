function map = DJK_hsv(m)
%   DJK_hsv : Hue-saturation-value color map
%   
%   071207: Adjusted from hsv(). In hsv() only the hue is varied in the 
%   colormap. In order to get more distinct colors, in this colormap the 
%   saturation and value (light) are also varied in big discrete steps

% DETERMINE SIZE OF COLORMAP
if nargin < 1, m = size(get(gcf,'colormap'),1); end

% HUE IS GRADIENT FROM 0 to 1
h = (0:m-1)'/max(m,1);

% SATURATION is either 0.4 0.7 or 1.0
temp = 0.4;
for i = [1:m]
    s(i,1) = temp;
    temp = temp + 0.3;
    if temp>1
        temp = 0.4;
    end
end

% VALUE is either 0.6 0.8 1.0 0.7 or 0.9
temp = 0.6;
for i = [1:m]
    v(i,1) = temp;
    temp = temp + 0.2;
    if temp>1
        temp = temp - 0.5;
    end
end

if isempty(h)
  map = [];
else
  map = hsv2rgb([h s v]);
end

