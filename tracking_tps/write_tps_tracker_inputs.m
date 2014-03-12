function [] = write_tps_tracker_inputs(integerImage,filename)

fid = fopen(filename,'w');

% Get the IDs of actual cells which occur in the image (0 is background)
% (to skip any integers that are not used in the image)
cellIds = setdiff(unique(integerImage),0)';

r = regionprops(integerImage,'Centroid','Orientation','MajorAxisLength');
cen = [r.Centroid];
cenx = cen(1:2:end);
ceny = cen(2:2:end);
o = [r.Orientation];
len = [r.MajorAxisLength];

% We don't have 3D data, only 2D, so set Z coordinate to 1 for all cells.
% cenz same length as other properties, to keep aligned with other properties
maxCellId = length(cenx);
cenz = ones([1,maxCellId]);

data = [cenx(cellIds); ceny(cellIds); cenz(cellIds); ...
        len(cellIds); o(cellIds); cellIds];
%data = [cenx(cellIds); ceny(cellIds); ...
%        len(cellIds); o(cellIds); cellIds];

fprintf(fid,'# x, y, z, length, angle, cell_id\n');
fprintf(fid,'%f %f %f %f %f %d\n',data(:,:));
%fprintf(fid,'# x, y, length, angle, cell_id\n');
%fprintf(fid,'%f %f %f %f %d\n',data(:,:));
fclose(fid);

return;
