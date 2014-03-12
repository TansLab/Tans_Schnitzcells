
function treeend=drawschnitzbaumST(p,schnitzcells,me,x,level,treeend);

% draws a tree from a schnitzcell structure, starting with "me" as the
% root.
% when you run this yourself, omit "x" and "level" (those parameters are 
% used by drawschnitzbaum when it calls itself recursively.


miny = 0;
minc = 0;

ccolor = [1 1 1];
ycolor = [1 1 1];

cedgecolor = 0.5*[0 1 0];
yedgecolor = 0.5*[1 0 0];

% maxy = 3/2;
% maxc = 211/2;

maxc = 0.9*max([schnitzcells.MC]);
maxy = 0.9*max([schnitzcells.MY]);
% maxy = 0.9*max([schnitzcells(me).MY]);
maxy = 0.5*max([schnitzcells.MY]);

if nargin == 3,
    x = 0;
    level = 1;
    treeend=[];
    figure;
end;

epsilon = 0.01;
epsilon = 0.02;

% first draw ourself
frames = schnitzcells(me).frames;
if schnitzcells(me).approved
    plot(x*ones(length(frames),1),frames,'.-');
else
    plot(x*ones(length(frames),1),frames,'r.-');
end
hold on;
c = schnitzcells(me).MC;
y = schnitzcells(me).MY;
c(c<minc) = minc;
y(y<miny) = miny;

c = (c - minc)./(maxc-minc);
y = (y - miny)./(maxy-miny);

c(c>1) = 1;
y(y>1) = 1;

for i = 1:length(frames),
    if ~isempty(c),
        if ~isnan(c(i)),
            h1 = plot(x+epsilon,frames(i), 's');
            set(h1,'markerfacecolor', (c(i)*ccolor));
            set(h1,'markeredgecolor', cedgecolor);
            set(h1,'markersize', 10);
        end;
    end;
    if ~isempty(c),
        if ~isnan(y(i)),
            h2 = plot(x-epsilon,frames(i), 's');
            set(h2,'markerfacecolor', (y(i)*ycolor));
            set(h2,'markeredgecolor', yedgecolor);
            set(h2,'markersize', 10);
        end;
    end;
end;
% text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).cellno'),'color','r','fontsize',7);
%  text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).ancestor'),'color','r','fontsize',7);
% now draw our daughters:
%keyboard
if schnitzcells(me).D>0,
    text(x+1/2^level, frames(end), num2str(schnitzcells(me).D),'FontSize',6); %%ADD SJT 27-12-05
    line([x x+1/2^level],[frames(end) frames(end)+1])
    treeend=drawschnitzbaumST(p,schnitzcells,schnitzcells(me).D,x+1/2^level,level+1,treeend);
else 
    treeend = [treeend me];
end;
if schnitzcells(me).E>0,
    text(x-1/2^level, frames(end), num2str(schnitzcells(me).E),'FontSize',6); %%ADD SJT 27-12-05
    line([x x-1/2^level],[frames(end) frames(end)+1])
    treeend=drawschnitzbaumST(p,schnitzcells,schnitzcells(me).E,x-1/2^level,level+1,treeend);
elseif schnitzcells(me).D>0
    treeend = [treeend me];
end;

%keyboard
if nargin == 2,
%     set(gca,'color','k');
    a = axis;
    a(1:2) = [-1 1];
    axis(a);
end;

title([p.movieDate '  ' p.movieName]);


