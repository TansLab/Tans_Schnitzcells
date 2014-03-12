function newim = drawline(oldim,pt1,pt2,color)
% function newim = drawline(oldim,pt1,pt2,color);

vec = [pt2(1) - pt1(1) , pt2(2)-pt1(2)];
D = sqrt(sum(vec.^2));
if (D==0)
  oldim(pt1(1),pt1(2)) = color;
else
  for d = 0:0.25:D,
    thispt = round(pt1 + vec*d/D);
    oldim(thispt(1), thispt(2)) = color;
  end
end
newim = oldim;
