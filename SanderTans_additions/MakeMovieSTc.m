function mov=MakeMovieST(p,frstart,frend,PhOffset,FyOffset)

%
%
%
%
% for fr = frames,
% 
%     ind = ind + 1;
% 
%     % load the seg'd image:
% 
%     clear Lc rect; [Lc,rect] = loadseg(p,fr,'Lc','rect');

clear mex;
fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
       'NextPlot','replace','Visible','off')
basepathPh = [p.imageDir,p.movieName,'-p-'];
% if exist(p.numphaseslices)
%     if p.numphaseslices == 1
%         basepathPh = [p.imageDir,p.movieName,'-p-'];
%     else
%         % assumes phase image 2 is the one to extract from!  >> saves both
%         % slices SJT
%         basepathPh = [p.imageDir,p.movieName,'-p-'];
%         basepathPh2 = [p.imageDir,p.movieName,'-p-2-'];
%     end
% end

% bgDir=[p.rootDir,'Background\bg2.tif'];
% BG=imread(bgDir);

basepathFy = [p.imageDir,p.movieName,'-c-'];

moviefile = [p.imageDir,'PhFc2.avi'];
moviefile
p.imageDir
mov = avifile(moviefile,'fps',4,'compression','none','quality',100);
% mov = avifile(moviefile,'fps',3);

for i=frstart:frend
%     imp = imread([basepathPh,str3(i),'.tif']);
    imp = imread([basepathPh,'1-',str3(i),'.tif']);
    imy = imread([basepathFy,str3(i),'.tif']);    
%     imy = (imy - FyOffset - BG); %% ONLY IF BACKGRND SUBTRACT!!!
    imy = (imy - FyOffset); %% ONLY IF BACKGRND SUBTRACT!!!
    i
    maxmax(imy)
    minmin(imy)
    imshow(makergb(imresize(imp,0.5)-PhOffset,imy));
    %imshow(makergb(imresize(imp,1)-PhOffset,imy));
%     keyboard
    F = getframe(gca);
    mov = addframe(mov,F);
end
mov = close(mov);
clear mex;