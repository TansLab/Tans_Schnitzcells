function treeend=gettreeendST(p,schnitzcells,startschnitz)
%
%
%
treeend=[];
for me=startschnitz:length(schnitzcells)
    if (schnitzcells(me).P>0 & schnitzcells(me).D==0 & schnitzcells(me).E==0),
        treeend = [treeend me];
    end;
end
