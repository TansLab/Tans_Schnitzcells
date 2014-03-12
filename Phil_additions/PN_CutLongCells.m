function cutterImage = PN_CutLongCells(thresholdCut,Image)
%outputs an binary image with positions to cut the skeleton
cutterImage = false(size(Image));
pixelsToSuppressXY = [];

dj0 = bwdist(~Image);
dj1 = imhmax(dj0,1);         
i1 = bwlabel(Image);

stats = regionprops(i1,'Area','BoundingBox','Image');
characSize = median([stats.Area]);
idx = find([stats.Area] > characSize*1.5); %suspicious cells

for ii = idx  %study the case of long cells individually
    s = stats(ii);
    xb = ceil(s.BoundingBox(1)); 
    yb = ceil(s.BoundingBox(2));
    lx = s.BoundingBox(3);
    ly = s.BoundingBox(4);
    s.Image(1,:) = 0; s.Image(:,1) = 0; s.Image(ly,:) = 0; s.Image(:,lx) = 0; 
    skelet = bwmorph(s.Image,'skel',inf); 
    skelet = bwmorph(skelet,'spur',20);   
    subImage = imcrop(dj1,[xb yb lx-1 ly-1]);
    subImage(~skelet) = inf;                
    [m xm ym] = MinCoordinates2(subImage);
    if m < thresholdCut
        pixelsToSuppressXY = [pixelsToSuppressXY ; xm+xb-1 ym+yb-1];
    end
end

if ~isempty(pixelsToSuppressXY)
    pixelsToSuppress = sub2ind(size(Image),pixelsToSuppressXY(:,2),pixelsToSuppressXY(:,1));
    cutterImage(pixelsToSuppress) = true;
end