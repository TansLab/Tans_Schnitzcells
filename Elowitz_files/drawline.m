function newim = drawline(oldim,pt1,pt2,color)
% function newim = drawline(oldim,pt1,pt2,color);
% draws line based, using # elements based on length of line

%%
vec = [pt2(1) - pt1(1) , pt2(2)-pt1(2)];
D = sqrt(sum(vec.^2));
if (D==0)
  oldim(pt1(1),pt1(2)) = color;
else        
  for d = 0:.25:D, 
    thispt = round(pt1 + vec*d/D);
    oldim(thispt(1), thispt(2)) = color;    
    %figure(4); hold on; plot(thispt(1), thispt(2),'.'); % debug feat
  end
end
newim = oldim;



