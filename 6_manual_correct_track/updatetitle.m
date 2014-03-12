function updatetitle;

global mousepos theimLc mainfig curschnitz dispframes curapproved;

numstr='';

xy = round(mousepos(1,1:2));
% [mfr, mx, my] = mouse2coords(xy, winsize, numcols);
if all(xy>0) & (xy(2)<size(theimLc,1)) & (xy(1)<size(theimLc,2)),
    cellno = theimLc(xy(2),xy(1));
    %     name = get(mainfig,'name');
    cellnumstr = '; cell # '; cellnumstrlen=length(cellnumstr);
    %     frnumstr = '; frame # '; frnumstrlen=length(frnumstr);
    %     f = findstr(name,cellnumstr);
    numstr = num2str(cellno);
    %     frstr = num2str(mfr);
    %     if isempty(f),
    %         name = cat(2,name,cellnumstr,numstr);
    %     else
    %         name = cat(2,name(1:f+cellnumstrlen),numstr);
    %     end;
    %     set(mainfig,'name',name);
end;
if curapproved,
    appchar = '*';
else
    appchar = '-';
end;

set(mainfig,'name',[appchar,' Schnitz: ',...
    num2str(curschnitz),', frames ',num2str(dispframes(1)),...
    '...',num2str(dispframes(end)),'; cell: ',numstr]);
pause(0.00001);
