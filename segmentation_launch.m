clear all
close all
vimages = [1:506];
%1 - cropping
p = DJK_initschnitz('pos9','2011-10-18','e.coli.AMOLF','rootDir','D:\colonies\', 'cropLeftTop', [1,1], 'cropRightBottom', [1200,1040]);
DJK_cropImages_2colors(p, vimages, [305,95], [1052,840], 'cropName', 'pos9crop');

p = DJK_initschnitz('pos9crop','2011-10-18','e.coli.AMOLF','rootDir','D:\colonies\', 'cropLeftTop', [305,95], 'cropRightBottom', [1052,840]);

%2 - segmentation
DJK_segLoopImproved_2colors(p, vimages, [1 2 3], [1 2 3], [275], 5);
DJK_copySegFilesImproved(p, vimages, [1 2 3], [1 2 3], [275], 5);

combinations = {[1 2 3],[1 2],[2 3]}

for ii=1:length(combinations)
    v = combinations{ii}
DJK_segLoopImproved_2colors(p,vimages, v, v, [275], 5);
end
DJK_copySegFilesImproved(p,vimages, [1 2 3], [1 2 3], [275], 5);
