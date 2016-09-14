
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

% Need negatives later on, so convert to signed int
my_image = int16(my_image);

% Show image
h1 = figure();
imshow(my_image,my_colormap);

% Select average from which to determine average background________________

% Select region from which to determine average
my_rect = round(getrect(h1))
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

m = figure();
hist(double(my_avg_area(:)),[double(min(my_avg_area(:))):1:double(max(my_avg_area(:)))]);

% Determine delta image____________________________________________________
% Delta image is value above mean for each pixel

my_delta_image = my_image;
my_delta_image = my_delta_image - my_mean;

h102 = figure();

%norm_my_delta_image = (my_delta_image - min(my_delta_image(:))) / (max(my_delta_image(:)) - min(my_delta_image(:)));

imshow(my_delta_image*-1,my_colormap)

% Determine what is dust and remove________________________________________

% Create black image of same size original
my_dust_image = zeros( size (my_image) );
% Mark area used for averaging
my_dust_image(my_ymin:my_ymax,my_xmin:my_xmax) = .5;
% Mark area to be left alone
%my_dust_image(my_notaltered_ymin:my_notaltered_ymax,my_notaltered_xmin:my_notaltered_xmax) = .75;

% Copy original image to to-be-dusted image
my_shiny_image = my_image;

% Identify dust and remove
for y = [1:size(my_image,1)]
    for x = [1:size(my_image,2)]
               
        % If outside region not to be altered (dusted)
        if ~(my_BW_not_altered(y,x))
        
            % Substract delta image 
            my_shiny_image(y,x) = my_image(y,x) - my_delta_image(y,x);
            
            % Detect dust and remove main body of dust
            %{
            % Remove dust
            if (my_image(i,j)<my_min_margin)
               % Mark place as dusty
               my_dust_image(i,j) = 1; 
               % Replace dusty value for avg 
               my_shiny_image(i,j)=my_mean;
            end
            %}
            
        end
        
    end
end

% Show dust image
h2 = figure();
imshow(my_dust_image);

% Save image_______________________________________________________________


% Convert back to uint
my_shiny_image = uint16(my_shiny_image);

% Show result
h3 = figure();
imshow(my_shiny_image,my_colormap);

% Get old image metadata
im_info = imfinfo([my_dir img_name]);

% Note that this produces 'cropped' image style metadata
im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
imwrite(my_shiny_image,[my_dir output_name], 'tif', 'Compression', 'none', 'Description', im_description);


%%
% Now also for other images________________________________________________

% configs
% offsets need to be determined for each image
xoffset = -26;
yoffset = -57;

% Load image_______________________________________________________________
my_dir = 'D:\LocalPlayground\'
img_name = 'pos1-p-3-101.tif'
output_name = ['D' img_name]
%img_location = 'D:\LocalPlayground\Vincent_van_Gogh_(1853-1890)_-_Wheat_Field_with_Crows_(1890).jpg'
[my_image,my_colormap] = imread([my_dir img_name]);

% apply offset - probably also possible by fn
offset_my_image = ones(size(my_image))*my_mean; % TODO wiped later, so useless

source_x1 = max(1,1-xoffset) 
source_y1 = max(1,1-yoffset)
source_x2 = min(size(my_image,2),size(my_image,2)-xoffset)
source_y2 = min(size(my_image,1),size(my_image,1)-yoffset)

target_x1 = max(1,1+xoffset) 
target_y1 = max(1,1+yoffset)
target_x2 = min(size(my_image,2),size(my_image,2)+xoffset)
target_y2 = min(size(my_image,1),size(my_image,1)+yoffset)

offset_my_image([target_y1:target_y2],[target_x1:target_x2]) = my_image([source_y1:source_y2],[source_x1:source_x2]);
my_image = offset_my_image;

% Need negatives later on, so convert to signed int
my_image = int16(my_image);

% Show image
h1 = figure();
imshow(my_image,my_colormap);
 
% Subtract average image __________________________________________________

% Copy original image to to-be-dusted image
my_shiny_image = my_image;

% boundary indices that fall outside comparable ranges of images
bx1 = abs(xoffset); by1 = abs(yoffset);
bx2=size(my_shiny_image,2)-abs(xoffset); by2 = size(my_shiny_image,1)-abs(yoffset);

% Identify dust and remove
for y = [1:size(my_image,1)]
    for x = [1:size(my_image,2)]
               
        % If outside region not to be altered (dusted)
        if ~(my_BW_not_altered(y,x))
        
            % Substract delta image 
            my_shiny_image(y,x) = my_image(y,x) - my_delta_image(y,x);
            
            % Detect dust and remove main body of dust
            %{
            % Remove dust
            if (my_image(i,j)<my_min_margin)
               % Mark place as dusty
               my_dust_image(i,j) = 1; 
               % Replace dusty value for avg 
               my_shiny_image(i,j)=my_mean;
            end
            %}
            
        end
        
        if x < bx1 | x > bx2 | y < by1 | y > by2
            my_shiny_image(y,x) = my_mean;
        end
        
    end
end

% print____________________________________________________________________

% Convert back to uint
my_shiny_image = uint16(my_shiny_image);

% Show result
h3 = figure();
imshow(my_shiny_image,my_colormap);

% Get old image metadata
im_info = imfinfo([my_dir img_name]);

% Note that this produces 'cropped' image style metadata
im_description = [im_info.ImageDescription 'DateTime: ' im_info.DateTime 'Software: ' im_info.Software];
imwrite(my_shiny_image,[my_dir output_name], 'tif', 'Compression', 'none', 'Description', im_description);



