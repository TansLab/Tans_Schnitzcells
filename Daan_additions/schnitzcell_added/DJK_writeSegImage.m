function DJK_writeSegImage(L,filename);
% DJK_writeSegImage is used to write a png image of segmentation image 

if islogical(L),
  L = double(L);
end;

% L2 has every non-background blob in range [2,256] and 
% sets background to one, corresp. to first entry in mymap
L2 = mod(L,255)+2;
L2(L==0) = 1;

% M is the maximum color table entry, at most 256 colors
M = min(max2(L)+2,256);
% create a color map
mymap = DJK_hsv(M); %DJK 071207
% explicitly set the colormap's first entry to black for background
mymap(1,:)=[0 0 0];
% get sequence of random integers in range [1,maxcolors-1]
[s,I] = sort(rand(M-1,1));  
% randomly reorder mymap color entries [2,maxcolors]
mymap(2:end,:) = mymap(I+1,:);

imwrite(L2,mymap,filename,'png')
