vimages = [77];

DJK_segLoopImproved_2colors(p,vimages, [1], [1], [275], 5);

DJK_copySegFilesImproved(p,vimages, [1], [1], [275], 5);
DJK_manualcheckseg(p,'manualRange',vimages,'override',0);
