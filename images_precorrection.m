clear all
close all

imageDir = 'X:\colonies\2011-10-07\pos3crop\images\'

D = dir([imageDir, ['*.tif'] ]);

for ii = 1:length(D)
    D(ii).name = [ D(ii).name(1:4) 'crop' D(ii).name(5:length(D(ii).name))];
end



%    
%     for i = [1:length(DphaseRange)], % go over each phase image of this frame
%         % read image
%         im_original = imread([imageDir DphaseRange(i).name]); 
%         % get image info
%         im_info = imfinfo([imageDir DphaseRange(i).name]);
%         % this image info will be added to crop
%         im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
%         % crop the image
%         im_crop = im_original(:,1:1200);
%         % write image data
%         im_crop_filename = [DphaseRange(i).name];
%         imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
%         disp(['Written: ' im_crop_filename]);
%         %disp([ num2str(leftTop(2)) ':' num2str(rightBottom(2)) '-' num2str(leftTop(1)) ':' num2str(rightBottom(1))]);
%     end
%     
%     DfluorYRange = dir([imageDir, ['pos3-y-*'] ]);
%     DfluorRRange = dir([imageDir, ['pos3-r-*'] ]);
%     for i = [1:length(DfluorYRange)],  % go over each fluor image of this frame
%         % read image
%         im_original = imread([imageDir DfluorYRange(i).name]); 
%         % get image info
%         im_info = imfinfo([imageDir DfluorYRange(i).name]);
%         % this image info will be added to crop
%         im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
%         % crop the image
%         im_crop = im_original(:,1:(1200/2));
%         % write image data
%         im_crop_filename = [DfluorYRange(i).name];
%         imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
%         disp(['Written: ' im_crop_filename]);
%         %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
%     end
%     for i = [1:length(DfluorRRange)],  % go over each fluor image of this frame
%         % read image
%         im_original = imread([imageDir DfluorRRange(i).name]); 
%         % get image info
%         im_info = imfinfo([imageDir DfluorRRange(i).name]);
%         % this image info will be added to crop
%         im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
%         % crop the image
%         im_crop = im_original(:,1:(1200/2));
%         % write image data
%         im_crop_filename = [DfluorRRange(i).name];
%         imwrite(im_crop, im_crop_filename, 'tif', 'Compression', 'none', 'Description', im_description);
%         disp(['Written: ' im_crop_filename]);
%         %disp([ num2str(leftTop(2)/2) ':' num2str(rightBottom(2)/2) '-' num2str(leftTop(1)/2) ':' num2str(rightBottom(1)/2)]);
%     end

clear all