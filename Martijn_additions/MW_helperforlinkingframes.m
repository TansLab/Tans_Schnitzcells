

FRAMEN=155;

filename1 = [p.segmentationDir,p.movieName,'seg',str3(FRAMEN),'.mat'];
load(filename1);
LcFrame1 = Lc;
filename2 = [p.segmentationDir,p.movieName,'seg',str3(FRAMEN+1),'.mat'];
load(filename2);
LcFrame2 = Lc;

h = figure(1); clf; hold on;

linklist=[];
while(any(LcFrame1(:)>0))
    disp('.');
    
    MW_showtwoframeslabeled(LcFrame1, LcFrame2,h)

    click1 = round(ginput(1))
    cellno1 = LcFrame1(click1(2),click1(1))
    click2 = round(ginput(1))
    cellno2 = LcFrame2(click2(2),click2(1))

    linklist = [linklist; cellno1, cellno2]
    
    LcFrame1(find(LcFrame1==cellno1)) = 0;
    LcFrame2(find(LcFrame2==cellno2)) = 0;
end





