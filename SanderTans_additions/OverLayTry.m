% Create the scene
% This example uses a random set of pixels to create a TrueColor image
scene = rand(100,100,3);

% Create the object image
% This example uses a blend of colors from left to right, converted to a TrueColor image
% Use repmat to replicate the pattern in the matrix
% Use the "jet" colormap to specify the color space
obj = ind2rgb(repmat(1:64,64,1),jet(64));

% Display the images created in subplots
hf2 = figure('units','normalized','position',[.2 .2 .6 .6]);
axi1 = subplot(2,3,1);
iscene = image(scene);
axis off
title('Scene')
axi2 = subplot(2,3,4);
iobj = image(obj);
axis off
title('Object image')

% Now replace pixels in the scene with the object image
result = scene;
% Define where you want to place the object image in the scene
rowshift = 20;
colshift = 0;
% Perform the actual indexing to replace the scene's pixels with the object
result((1:size(obj,1))+rowshift, (1:size(obj,2))+colshift, :) = obj;
% Display the results
ax3 = subplot(2,3,[2:3, 5:6]);
iresult = image(result);
axis off
hold on
title(sprintf('Using indexing to overlay images:\nresult is one image object'))% Let's add another image to the scene. 
% First define the polygon region of interest
t = 0:.1:2*pi;
[objheight,objwidth] = size(obj(:,:,1));
x = round(20*cos(t)+round(objwidth/2));
y = round(20*sin(t)+round(objheight/2));
% Since the mask is 2D, repeat it through the 3D RGB matrix
% while extracting the polygon mask
objCircleMask = repmat(roipoly(obj,x,y),[1,1,3]);
% Define where we want to place the circle object in the scene
rowshift = 0;
colshift = 30;
result2 = result;
% Perform the actual indexing. First calculate the replacement pixels
% by multiplying the scene by the negative of the mask, and then adding 
% the object multiplied by the mask. Multiplying the mask chooses all 
% pixels corresponding to a 1 in the mask, and removes color from the 
% pixels corresponding to a 0.
replacementPixels = ( result((1:size(obj,1))+rowshift, (1:size(obj,2))+colshift, :) ...
.* double(~objCircleMask) ) + ( obj .* double(objCircleMask) );
% Now replace the pixels with the ones we just calculated
result2((1:size(obj,1))+rowshift, (1:size(obj,2))+colshift, :) = replacementPixels;
% Display the results. Notice that we placed the circle object in the scene
% that resulted from the previous object placement.
axes(ax3);
delete(iresult)

iresult2 = image(result2);
