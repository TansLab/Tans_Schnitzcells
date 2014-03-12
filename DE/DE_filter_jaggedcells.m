function DE_filter_jaggedcells(B_edgeImage2,howstrict)
%% 08 july 2013
% Some experiments result with separated cells. The empty space bewteen
% them very often is recognized as a cell, and it consumes quite some time
% to remove them manually.
% This routine is an attempt to throw a check on a region, taking the
% ratio between R=Perimeter/Area. Jageed those regions will have unusually 
% high R; thus we can filter them out by simple comparison to mean(R).
allcells=0;

zz=B_edgeImage2;
xx=imfill(zz,'holes');
cc=regionprops(xx,'Area','Perimeter','Centroid');
c1=[cc.Area]';
c2=[cc.Perimeter]';
c3=[];
for i=1:length(cc)
    c3(i,:)=cc(i).Centroid;
end

R=c2./c1;
R_avg=mean(R);
[ind,~]=find(R>=howstrict*R_avg);

figure;
hold on
imagesc(xx)
if allcells==1
    for i=1:length(cc)
        text(c3(i,1),c3(i,2),num2str(i))
    end
else
    for i=ind
        text(c3(i,1),c3(i,2),num2str(i))
    end
end

