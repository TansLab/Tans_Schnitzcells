function xpand_segimage(p,imnamenum,amount);
% function X_expand_segimage(imnamenum,amount);
% Incase the segmented image cuts off some cells. 
% Resaves the images of frame number imnamenum
% with bigger size expanded by amount.

% NOTES on variable name changes
%   tifdir    -> p.imageDir
%   imprefix  -> p.movieName
%   regsize   -> p.regsize
%   outprefix -> p.outprefix

mnamePhase2   = [p.movieName,'-p-2-*.tif'];
mnamePhaseAll = [p.movieName,'-p-*.tif'];

Dphase2   = dir([p.imageDir, mnamePhase2]);
DphaseAll = dir([p.imageDir, mnamePhaseAll]);

% This code is not robust: problem if there are 2 or 3 phases and you set 
% numphaseslices = 1, D is too big.  Smarter to figure out which D is 
% smaller, and use the smaller one (BTW, same code's in segmoviephase)
if p.numphaseslices == 1
  D = DphaseAll;
else
  D = Dphase2;
end

numpos= findstr(D(1).name, '.tif')-3;
% JCR: This code not quite ideal- imnamenum is index into directory, 
% so we're not able to set manualrange the way we want- by frame number; 
% FIXME: change this code to work with an image number, not directory index num
name= D(imnamenum).name
% name(numpos-2) = '2';
% name(numpos:numpos+2) = str3(imnamenum);
mynum= name(numpos:numpos+2);
%mynum = str3(imnamenum);
pname = [p.imageDir,name];
phim = imread(pname);
Lname= [p.outprefix,mynum];   
cname= [p.imageDir,p.movieName,'-c-',mynum,'.tif'];
yname= [p.imageDir,p.movieName,'-y-',mynum,'.tif'];
gname= [p.imageDir,p.movieName,'-g-',mynum,'.tif'];

load([p.segmentationDir,Lname])
newrect = [ max(1+p.regsize,rect(1)-amount), ...
	    max(1+p.regsize,rect(2)-amount), ...
	    min(size(phim,1)-p.regsize,rect(3)+amount), ...
	    min(size(phim,2)-p.regsize,rect(4)+amount)      ];
newLNsub = zeros(size(phim));
newLNsub(rect(1):rect(3), rect(2):rect(4)) = LNsub;
rect = newrect;
LNsub = newLNsub(rect(1):rect(3), rect(2):rect(4));
phsub = phim(rect(1):rect(3), rect(2):rect(4));
savelist=['''phsub'',''LNsub'',''rect'''];
if exist(cname)==2
    disp('found CFP image');
    [creg, cshift, cback]= quickreg(LNsub, cname, rect, p.regsize, p.fullsize);
    [exptcstr, gainc, exptc]= imsettings(cname);
    savelist=[savelist,',''creg'',''cshift'',''exptc'',''gainc'',''cback'''];
end 
if exist(yname)==2
    disp('found YFP image');

    [yreg, yshift, yback]= quickreg(LNsub, yname, rect, p.regsize, p.fullsize);
    [exptystr, gainy, expty]= imsettings(yname);
    savelist=[savelist,',''yreg'',''yshift'',''expty'',''gainy'',''yback'''];
end 
if exist(gname)==2
    disp('found GFP image');
    [greg, gshift, gback]= quickreg(LNsub, gname, rect, p.regsize, p.fullsize);
    [exptystr, gaing, exptg]= imsettings(gname);
    savelist=[savelist,',''greg'',''gshift'',''exptg'',''gaing'',''gback'''];
end 

eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
disp(['saved file ',p.segmentationDir,Lname]);