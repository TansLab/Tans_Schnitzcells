function [xout,yout,xavgout,yavgout] = multiplotST(schnitzcells, xfield, yfield, whichones, style)
%
%
%
%
[x,y,xav,yav] = plotschnitzme(schnitzcells, xfield, yfield, whichones, style);

xm=cell2mat(x');
ym=cell2mat(y');
xavm=cell2mat(xav');
yavm=cell2mat(yav');

figure;
plot(xm',ym');


xout=xm;
yout=ym;
xavgout=xav;
yavgout=yav;

hold on;
derivon=0;
normon=0;
fs=findstr(yfield,'_');
for i=fs(end:-1:1);
    yfield(i+1:end+1)=yfield(i:end);yfield(i)='\';end
fs=findstr(xfield,'_');
for i=fs(end:-1:1);
    xfield(i+1:end+1)=xfield(i:end);xfield(i)='\';end
if derivon,
    ylabel(['derivative: \partial(',yfield,')/\partial(',derivfield,')']);
elseif normon,
    ylabel(['normalized ',yfield]);
else
    ylabel(yfield);
end;
xlabel(xfield);
hold off;


end