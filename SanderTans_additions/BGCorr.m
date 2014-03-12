%
%
% to substract Background images
%
%
p = initschnitz('Pos2','2005-15-11','e.coli','rootDir','D:\Projects\Schnitzcell');
bg1Dir=[p.rootDir,'Background\bg1.tif'];
bg2Dir=[p.rootDir,'Background\bg2.tif'];
bg3Dir=[p.rootDir,'Background\bg3.tif'];
bg4Dir=[p.rootDir,'Background\bg4.tif'];
bg5Dir=[p.rootDir,'Background\bg5.tif'];
BG1=imread(bg1Dir);
BG2=imread(bg2Dir);
BG=(BG1+BG2)/2;
BG3=imread(bg3Dir);
BG=(BG+BG2)/2;
BG4=imread(bg4Dir);
BG=(BG+BG2)/2;
BG5=imread(bg5Dir);
BG=(BG+BG2)/2;

im=imread([p.imageDir,'Pos2-y-001.tif']);

im2=(im+200);
im3=im2-BG2;
%imwrite([p.rootDir,'Background\bg1.tif']);