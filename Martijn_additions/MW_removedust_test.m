
% Configuration____________________________________________________________
% Used
min_multiplier = 0.9
% Not used
my_prctile = .01; % Percentile of avg img to be used for removal of dust

% Load image_______________________________________________________________
my_dir = 'D:\LocalPlayground\'
img_name = 'pos1-p-1-002.tif'
output_name = ['D' img_name]
%img_location = 'D:\LocalPlayground\Vincent_van_Gogh_(1853-1890)_-_Wheat_Field_with_Crows_(1890).jpg'
[my_image,my_colormap] = imread([my_dir img_name]);

% Show image
h = figure(1);
imshow(my_image,my_colormap);

% Select average from which to determine average background________________

% Select region from which to determine average
my_rect = round(getrect(h))
% Rename vars for clearity - note counter-intuitive x/y 
my_xmin = my_rect(1); my_ymin = my_rect(2); 
my_xmax = my_rect(1)+my_rect(3); my_ymax = my_rect(2)+my_rect(4);

% Select region which shouldn't be altered
my_BW_not_altered = roipoly;

% Get area
my_avg_area = my_image(my_ymin:my_ymax,my_xmin:my_xmax);
% Get average
my_mean = mean(my_avg_area(:));
% Get maximum
my_min = min(my_avg_area(:));
my_min_margin = my_min * min_multiplier;
my_min_prctile = prctile(my_avg_area(:),my_prctile);

m = figure(101);
hist(double(my_avg_area(:)),[double(min(my_avg_area(:))):1:double(max(my_avg_area(:)))]);


% Determine what is dust and remove________________________________________

% Create black image of same size original
my_dust_image = zeros( size (my_image) );
% Mark area used for averaging
my_dust_image(my_ymin:my_ymax,my_xmin:my_xmax) = .5;
% Mark area to be left alone
my_dust_image(my_notaltered_ymin:my_notaltered_ymax,my_notaltered_xmin:my_notaltered_xmax) = .75;

% Copy original image to to-be-dusted image
my_shiny_image = my_image;

% Identify dust and remove
for i = [1:size(my_image,1)]
    for j = [1:size(my_image,2)]
               
        % If outside region not to be altered (dusted)
        if ~(my_BW_not_altered(i,j))
        
            % Remove dust
            if (my_image(i,j)<my_min_margin)
               % Mark place as dusty
               my_dust_image(i,j) = 1; 
               % Replace dusty value for avg 
               my_shiny_image(i,j)=my_mean;
            end
            
        end
        
    end
end

% Show dust image
h = figure(2);
imshow(my_dust_image);

% Remove dust______________________________________________________________





% Show result
h = figure(3);
imshow(my_shiny_image,my_colormap);

% Get old image metadata
im_info = imfinfo([my_dir img_name]);

% Note that this produces 'cropped' image style metadata
im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
imwrite(my_shiny_image,[my_dir output_name], 'tif', 'Compression', 'none', 'Description', im_description);


 








