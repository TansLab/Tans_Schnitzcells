function outim = imshowlabel(L,bwscreen);
% IMSHOWLABEL is used to display an integer image 

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
mymap = hsv(M);
% explicitly set the colormap's first entry to black for background
mymap(1,:)=[0 0 0];
% get sequence of random integers in range [1,maxcolors-1]
[s,I] = sort(rand(M-1,1));  
% randomly reorder mymap color entries [2,maxcolors]
mymap(2:end,:) = mymap(I+1,:);

if nargin>=2,
	rgb = 0.5 * ind2rgb(L2,mymap);
	bwscreen = double(bwscreen);
	bwscreen = 0.5 * bwscreen / maxmax(bwscreen);
	rgb(:,:,1) = rgb(:,:,1) + bwscreen;
	rgb(:,:,2) = rgb(:,:,2) + bwscreen;
	rgb(:,:,3) = rgb(:,:,3) + bwscreen;
	imshow(rgb);
	outim = rgb;
else
	imshow(L2, mymap);	
    outim = L2;
end;
