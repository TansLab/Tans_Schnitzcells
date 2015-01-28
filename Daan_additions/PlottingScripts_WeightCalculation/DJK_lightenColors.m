% DJK_lightenColors returns lighter versions of given RGB colors
% * optional to add lightning as value from 1 to 100 (default 50)
%
% Code written by Daan Kiviet

function colors_out = DJK_lightenColors(colors, lightning)

if nargin == 1
    lightning = 50;
end

for i=1:length(colors)
  reds    = linspace(colors{i}(1), 1);
  greens  = linspace(colors{i}(2), 1);
  blues   = linspace(colors{i}(3), 1);
  colors_out{i} = [reds(lightning) greens(lightning) blues(lightning)];
end  
