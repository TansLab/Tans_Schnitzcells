

DJK_segLoopImproved_2colors(p, [121 391], [1 2], [1 2], [275], 5);
DJK_copySegFilesImproved(p, [121 391], [1 2], [1 2], [275], 5);

%3 - check / correct segmentation
DJK_manualcheckseg(p,'manualRange',[121 391],'override',0);


