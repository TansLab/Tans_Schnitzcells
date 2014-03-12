DJK_analyzeSeg(p,'manualRange',[1:720]);
DJK_tracker_djk(p,'manualRange', [1:720]);
problems = DJK_analyzeTracking(p,'manualRange', [1:720], 'pixelsMoveDef', 10, 'pixelsLenDef', [-4 6]);
%DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1,'problemCells',problems);
%if no problem anymore
DJK_makeMovie (p, 'tree', 'schAll', 'stabilize', 1);


