function DJK_drawschnitzbaum(p,schnitzcells,me,x,level);

% draws a tree from a schnitzcell structure, starting with "me" as the
% root.
% when you run this yourself, omit "x" and "level" (those parameters are 
% used by drawschnitzbaum when it calls itself recursively.


% IF THIS IS NOT A RECURSIVE RUN, SOME SETTINGS
if nargin == 3,
    x = 0;
    level = 1;
    h = figure;
    %maximize(h); % MW edit
end;


% COLOR SETTINGS
miny = 0;
maxy = 0.9*max([schnitzcells.Y6_mean_all]); %0.02
% yedgecolor = [0 0 0];


% OFFSET OF BOXES AND TEXT
offsetx_box = 0.00;
offsetx_text = 0.01;
offsety_text = 5;


% GET DATA
mins = schnitzcells(me).time;


% DRAW OURSELF - LINES
plot(x*ones(length(mins),1),mins,'k.-');
hold on;
text(x-offsetx_text, mins(1)-offsety_text, num2str(me),'FontSize',6); %% DJK 071204


% DRAW OUR DAUGHTERS AND CONNECTORS
if schnitzcells(me).D>0,
    mins_daughter = schnitzcells(schnitzcells(me).D).time;
    l2 = line([x x+1/2^level],[mins(end) mins_daughter(1)]);
    set(l2,'Color',[0 0 0]);
    DJK_drawschnitzbaum(p,schnitzcells,schnitzcells(me).D,x+1/2^level,level+1);
end;
if schnitzcells(me).E>0,
    mins_daughter = schnitzcells(schnitzcells(me).E).time;
    l2 = line([x x-1/2^level],[mins(end) mins_daughter(1)]);
    set(l2,'Color',[0 0 0]);
    DJK_drawschnitzbaum(p,schnitzcells,schnitzcells(me).E,x-1/2^level,level+1);
end;



% DRAW OURSELF - COLORS BOXES
y = schnitzcells(me).Y6_mean_all;
y(y<miny) = miny;
y = (y - miny)./(maxy-miny);
y(y>1) = 1; %y is now a value between 0 and 1, depending on MY and set miny and maxy
for i = 1:length(mins),
    if ~isnan(y(i)),
        h2 = plot(x-offsetx_box,mins(i), 's');
        set(h2,'markerfacecolor', DJK_getColor(y(i)));
        set(h2,'markeredgecolor', DJK_getColor(y(i))); %yedgecolor
        set(h2,'markersize', 10);
    end;
end;


% ADD AXIS AND COLOR BAR AND ASK TO SAVE
if nargin == 3,
    xlim([-1 1]);
    ylabel('Time (mins)');
    limit_yaxis = ylim;
    ylim([(mins(1)-10) limit_yaxis(2)]);
   
    cBar = colorbar('location','southoutside');
    %set(cBar,'YLim',[0 1]);
    maxMY = max([schnitzcells.Y6_mean_all]);
    set(cBar,'XTickLabel', {miny, 'max', maxy}); %((maxy-miny)/2)
    maxMY_index = round(1+64*maxMY/(maxy-miny));
    set(cBar,'XTick', [0 65 maxMY_index]);
    %get(cBar)
    
    %cBar_pos=get(cBar,'Position');
    %cBar_pos(3)=0.03;
    %set(cBar,'Position',cBar_pos)
    
    % ASK TO SAVE THE FIGURE
    saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
    pause(0.2);
    if (upper(saveFigInput(1))=='Y')
        saveas(gcf,[p.analysisDir p.movieName '-tree-' num2str(me) '.fig']);
        saveSameSize(gcf,'file',[p.analysisDir p.movieName '-tree-' num2str(me) '.png'], 'format', 'png');
        if strcmp(saveFigInput,'Yes and Close')
            close(gcf);
            pause(0.2);
        end
    end
end;

end

function color = DJK_getColor(x);
    if (x>1) 
        index=64; 
    elseif (x<=0)
        index=1;
    else
        index = ceil(64*x);
    end
    cMap = colormap('jet');
    color = cMap( index , :);
end