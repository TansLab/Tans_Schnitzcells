function imageReady = NW_prepareImagesRich(ph3,slices)
% TODO!!! (if possible at all)

imageReadyTemp = mean( ph3(:, :, slices), 3);

%scale image, such that 25th highest value is 10000 and the 25th lowest value is 0
x = double(imageReadyTemp);
s = sort(x(:));
small = s(25);
big = s(end-25);
rescaled = (x - small)/(big - small);
rescaled(rescaled<0) = 0;
imageReady = uint16(10000*rescaled);
