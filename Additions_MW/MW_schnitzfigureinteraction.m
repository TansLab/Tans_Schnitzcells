
function MW_schnitzfigureinteraction(src,callbackdata)
% note redundancy with MW_schnitzfigurewaitforclick
% is this function still used?

global pos currentFrameStr Limage ourfig pp phfig;

% position in figure
pos=max(1,round(get(gca,'CurrentPoint')));

% obtain intensity from image
if (pos(1,2)>0 && pos(1,2)<size(Limage,1) && pos(1,1)>0 && pos(1,1)<size(Limage,2));
    curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));
else
    curr_val = '-1';
end;

% set name in the figure
set(ourfig,'name',...
    ['Frame ' currentFrameStr ', Pos: ', num2str(pos(1,2)), ' , ',num2str(pos(1,1)),'  Val: ',curr_val]);
    
end