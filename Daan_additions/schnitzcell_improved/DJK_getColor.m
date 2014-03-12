function color = DJK_getColor(x);
if (x>1) 
    index=64; 
elseif (x<=0)
    index=1;
else
    index = ceil(64*x);
end
cMap = colormap('jet');
color = cMap( index , :);