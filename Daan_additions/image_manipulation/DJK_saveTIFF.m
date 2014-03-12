function DJK_saveTIFF(image, imageFileName) 
% DJK_saveTIFF saves an image to harddisk
%
% Required input:
%   * image: image in <10x10 uint16> format 
%   * imageFileName: image filename string without extension
%     ('D:\Analysis_Microscope_Data\test')
%

imwrite(image, [imageFileName '.tif'], 'TIFF', 'Compression', 'none');
disp(['saved: ' imageFileName '.tif']);
