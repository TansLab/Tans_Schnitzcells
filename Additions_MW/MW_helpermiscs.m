
% find a schnitz in a certain segmentation
cellLabelMYNEEDLE = 67;
%cellLabelMYNEEDLE = 118;
Lselect=Lc; Lselect(Lc>0)=1; Lselect(Lc==cellLabelMYNEEDLE)=2;imshow(Lselect,[])