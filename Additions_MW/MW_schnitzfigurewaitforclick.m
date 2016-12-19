
function MW_schnitzfigurewaitforclick(src,callbackdata)

error('This function has been marked for deletion, use MW_schnitzfigureinteraction instead')

%{
['global pos Limage ourfig res pp phfig;pos=max(1,round((1/res)*get(gca,''CurrentPoint'')));',...
                'if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2));',...
                'curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));else;curr_val=''-1'';end;',...
                'set(ourfig,''name'',[''Pos: '',num2str(pos(1,2)),'' , '',num2str(pos(1,1)),',...
                '''  Val: '',curr_val]);']
%}

global pos Limage ourfig pp phfig;


pos=max(1,round(get(gca,'CurrentPoint')));

if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2))
    curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));
else
    curr_val='-1';
end;

set(ourfig,'name',['Pos: ',num2str(pos(1,2)),' , ',num2str(pos(1,1)),'  Val: ',curr_val]);
            
            
end