function  [Lout,OKorNot,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,finetuneimage,gotoframenum] = ...
    manual_kant(p,Lin,leftend,upend,phin,finetuneimage,rect,min_size,imnum,expandvalue);

% iptsetpref('imshowtruesize','manual');
iptsetpref('imshowborder','tight');
backwards = 0;
global pos Limage ourfig res pp phfig

pp=0;pps=0;
zd3=25;
pb=[3,3,3,4;4,4,3,4;1,2,1,2;1,2,2,1];
bb=0;bbs=0;

newrect=[0,0,0,0];
crop_pop=0;
dontsave = 0;
addtolist=0;
savetemp=0;
gotoframenum=0;

Limage=Lin;
Lout=Lin;
OKorNot=0;
done=0;
quit_now=0;
pos=[1 1];

while ~done
clear j*  %j1=0;j2=0;
figure(ourfig)
clf reset;
imshowlabel(imresize(Lout,res));
pos12=get(ourfig,'position');

set(ourfig,'position',[leftend,upend-pos12(4),pos12(3),pos12(4)]);

set(ourfig,'name',['Pos: ',num2str(pos(1,2)),' , ',num2str(pos(1,1)),...
                         '  Val: ',num2str(double(Limage(pos(1,2),pos(1,1))))]);
figure(ourfig)

set(ourfig,'WindowButtonMotionFcn',['global pos Limage ourfig res pp phfig;pos=max(1,round((1/res)*get(gca,''CurrentPoint'')));',...
   'if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2));',...
   'curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));else;curr_val=''-1'';end;',...
   'set(ourfig,''name'',[''Pos: '',num2str(pos(1,2)),'' , '',num2str(pos(1,1)),',...
   '''  Val: '',curr_val]);']);
%   '''  Val: '',curr_val]);pp=showbox(phfig,pos,pp,size(Limage),res);figure(ourfig);']);

realpress = 0;
while ~realpress,
    ct=waitforbuttonpress;
    cc=get(ourfig,'currentcharacter');
    if (ct==1) & isempty(cc),
        realpress=0;
    else
        realpress=1;
    end;
end;

set(ourfig,'WindowButtonMotionFcn','');
%if ct cc=get(1,'currentcharacter');else cc=get(1,'selectiontype');end
if ct 
    % cc=get(ourfig,'currentcharacter');
    if cc==' '
        OKorNot=1;
        done=1;
    elseif cc=='q'
        Lout=Lin;
        OKorNot=0;
        quit_now=1;
        done=1;
    elseif cc=='p'
        if pp(1) delete(pp(1));delete(pp(2));delete(pp(3));delete(pp(4));end;pp=0;
        if pps
            pps=0;
        else
            figure(phfig);
            zoomrect=[max([1,pos(1,2)-zd3]),min([size(Lout,1),pos(1,2)+zd3]),...
                    max([1,pos(1,1)-zd3]),min([size(Lout,2),pos(1,1)+zd3])]*res;
            for pc=1:4
                pp(pc)=line([zoomrect(pb(1,pc)),zoomrect(pb(2,pc))],[zoomrect(pb(3,pc)),zoomrect(pb(4,pc))]);
                set(pp(pc),'color','w');
                set(pp(pc),'LineWidth',2);
            end
            figure(ourfig);
            pps=1;
        end
    elseif cc=='b'
        if bb delete(bb);delete(bbp);end;bb=0;
        if bbs
            bbs=0;
        else
            cutx=round(pos(1,2));
            cuty=round(pos(1,1));
            if Lout(cutx,cuty)
                [fx,fy] = find(Lout==Lout(cutx,cuty));
                figure(phfig);
                bb=plot(fy*res,fx*res,'w.');
                set(bb,'markersize',2);
                bbp=plot(cuty*res,cutx*res,'r.');
                set(bbp,'markersize',24);
                figure(ourfig);
                bbs=1;
            end
        end
    elseif cc=='a'
        zoomrect=[max([1,pos(1,2)-50]),min([size(Lout,1),pos(1,2)+49]),...
                  max([1,pos(1,1)-50]),min([size(Lout,2),pos(1,1)+49])];
        Lzoom=Lout(zoomrect(1):zoomrect(2),zoomrect(3):zoomrect(4));
        Phzoom=phin(zoomrect(1):zoomrect(2),zoomrect(3):zoomrect(4));
        addfig=figure;
        
        LZedge = zeros(size(Lzoom));
        for ie = 1:max2(Lzoom);        
            LZedge = LZedge | bwperim(Lzoom==ie);
        end;
        LZedge=double(+LZedge);
        imshow(makergb(+imresize(LZedge,5),imresize(Phzoom(:,:,1),5)));
        figure(addfig); % needed in case current figure changes (in Windows)
        subaddcell=imresize(roipoly,1/5);%(phin);
        if max2(Lzoom(subaddcell>0))>0
            disp('overlaps existing cell; ignored.');
        else
            Lzoom(subaddcell)=max2(Lout)+1;
            Lout(zoomrect(1):zoomrect(2),zoomrect(3):zoomrect(4))=Lzoom;
        end
        close(addfig)
        figure(ourfig)
        done=0;
    elseif cc=='x'
        subcolroi=imresize(~roipoly,1/res);
        if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
            subcolroi2=zeros(size(Lout));
            subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
            subcolroi=subcolroi2;
        end
                Lout=(double(Lout).*double(subcolroi));
    %           Lout=renumberimage(Lout);
        done=0;
    elseif cc=='k'
        subcolroi=imresize(~roipoly,1/res);
        if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
            subcolroi2=zeros(size(Lout));
            subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
            subcolroi=~subcolroi2;
        end
                Lout=(double(Lout).*double(subcolroi));
    %           Lout=renumberimage(Lout);
        done=0;

    elseif cc=='t'

        subcolroi=imresize(~roipoly,1/res);
        if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
            subcolroi2=zeros(size(Lout));
            subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
            subcolroi=subcolroi2;
        end
                [xe,ye]=find(edge(subcolroi)==1);
                for qe=1:length(xe)
            if Lout(xe(qe),ye(qe))~=0
                Lout(Lout==Lout(xe(qe),ye(qe)))=0;
            end
                end
                Lout=(double(Lout).*double(subcolroi));
%               Lout=renumberimage(Lout);
        done=0;

    elseif cc == 'c'
        crop_pop=1;
        done=1;
                extra= 15; % <- extra number of pixels on either side of highlighted region
                LNfull=zeros(p.fullsize(1),p.fullsize(2));
                LNfull(rect(1):rect(3), rect(2):rect(4))=Lout;
                [fx,fy]= find(LNfull);
                xmin= max(min(fx) - extra, 1);
                xmax= min(max(fx) + extra, size(LNfull,1));
                ymin= max(min(fy) - extra, 1);
                ymax= min(max(fy) + extra, size(LNfull,2));
                newrect= [xmin ymin xmax ymax];
    elseif cc == 's'
        savetemp=1;
        done=1;
    elseif cc == 'e'
        disp('"expand image" currently not in use.');
%         expand_segimage(p,imnum,expandvalue);
%         OKorNot=1;
%         done=1;
%         dontsave=1;
%         backwards = 2;
    elseif cc == 'f'
        finetuneimage = 1;
    elseif cc == 'g'
        gotoframenum = input('goto loopindex = ');
        OKorNot=1;
        done=1;
        dontsave=1;
    elseif cc == 'w'
        savetemp=2;
        done=1;
    elseif cc == '.'
        OKorNot=1;
        done=1;
        dontsave=1;
   elseif cc == ','
        OKorNot=1;
        done=1;
        dontsave=1;
        backwards = 1;
    elseif cc == 'R'
        Lout=renumberimage(Lout);
    elseif cc == 'r'
        % renumber this cell
        cutx=round(pos(1,2));
        cuty=round(pos(1,1));
        chosencolor=Lout(cutx,cuty);
        cell=zeros(size(Lout));
        cell(Lout==Lout(cutx,cuty))=Lout(Lout==chosencolor);
        [fx,fy] = find(cell);
        xmin = max(min(fx)-5,1);
        xmax = min(max(fx)+5,size(cell,1));
        ymin = max(min(fy)-5,1);
        ymax = min(max(fy)+5,size(cell,2));
        subcell = cell(xmin:xmax, ymin:ymax);
        cell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4);    
        Lout(Lout==Lout(cutx,cuty))=0;
        Lout(cell==1) = chosencolor;
        for k = 2:max2(cell),
            Lout(cell==k) = max2(Lout)+k-1;
        end;
    elseif cc == 'l'
        addtolist=1;
        OKorNot=1;
        done=1;
        dontsave=1;
    elseif cc == 'o'
        % obliterate all but this cell
        cutx=round(pos(1,2));
        cuty=round(pos(1,1));
        chosencolor=Lout(cutx,cuty);
        if (chosencolor > 0)
          Lout = (Lout==chosencolor);
        end
    elseif cc==27
        Lout=Lin;OKorNot=0;
        done=1;
    elseif cc =='',
        disp('you typed shift');

    end
else 
    cz=get(ourfig,'selectiontype');
    % normal = left mouse button
        if cz(1)=='n' 
        % join cells
        pos1=pos;
        j1=Lout(round(pos(1,2)),round(pos(1,1)));
        figure(ourfig);
        set(ourfig,'WindowButtonMotionFcn',['global pos Limage ourfig res pp phfig;pos=max(1,round((1/res)*get(gca,''CurrentPoint'')));',...
            'if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2));',...
            'curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));else;curr_val=''-1'';end;',...
            'set(ourfig,''name'',[''Pos: '',num2str(pos(1,2)),'' , '',num2str(pos(1,1)),',...
            '''  Val: '',curr_val]);']);
        ct=waitforbuttonpress;
        set(ourfig,'WindowButtonMotionFcn','');
        if ct
            j2=j1;
        else
            j2=Lout(round(pos(1,2)),round(pos(1,1)));
        end
        if j2<j1
            j1old=j1;
            j1=j2;
            j2=j1old;
        end
        Lout(Lout==j2)=j1;
        if ((pos1(1,2)-pos(1,2))^2+(pos1(1,1)-pos(1,1))^2)
            Lout=drawline(Lout,[pos1(1,2),pos1(1,1)],[pos(1,2),pos(1,1)],j1);
            Lout=drawline(Lout,[pos1(1,2)+1,pos1(1,1)],[pos(1,2)+1,pos(1,1)],j1);
            Lout=drawline(Lout,[pos1(1,2),pos1(1,1)+1],[pos(1,2),pos(1,1)+1],j1);
        end
         % extend = shift+left button (or both?)
         elseif (cz(1)=='e' & Lout(round(pos(1,2)),round(pos(1,1))))
        % erase cell
        Lout(Lout==Lout(round(pos(1,2)),round(pos(1,1))))=0;
         % alternate) = ctrl+left button or right mouse button
         elseif (cz(1)=='a' & Lout(round(pos(1,2)),round(pos(1,1))))
        % cut at this location
        cutx=round(pos(1,2));
        cuty=round(pos(1,1));
        chosencolor=Lout(cutx,cuty);
        cell=zeros(size(Lout));
        cell(Lout==Lout(cutx,cuty))=Lout(Lout==chosencolor);
        [fx,fy] = find(cell);
        xmin = max(min(fx)-5,1);
        xmax = min(max(fx)+5,size(cell,1));
        ymin = max(min(fy)-5,1);
        ymax = min(max(fy)+5,size(cell,2));
        subcell = cell(xmin:xmax, ymin:ymax);

        perim=bwperim(imdilate(subcell,strel('disk',1)));
%        figure(1);imshowlabel(perim);pause;close(1);
        perims=zeros(size(perim));
        radp=1;
        while max2(perims)<2 & radp<41
            pxmin = max(cutx-xmin+1-radp,1);
            pxmax = min(cutx-xmin+1+radp,size(perims,1));
            pymin = max(cuty-ymin+1-radp,1);
            pymax = min(cuty-ymin+1+radp,size(perims,2));
            perims(pxmin:pxmax,pymin:pymax)=bwlabel(perim(pxmin:pxmax,pymin:pymax));
            radp=radp+1;
        end
        if max2(perims)>1
            kim=zeros(size(subcell));
            kim(cutx-xmin+1,cuty-ymin+1)=1;
            kim1=kim;
            while ~any(any(kim1 & perims))
                kim1=imdilate(kim1,strel('disk',1));
            end
            [cut1x,cut1y]=find(kim1 & perims);
            color1=perims(cut1x(1),cut1y(1));
            perims(perims==color1)=0;
            kim2=kim;
            while ~any(any(kim2 & perims))
                kim2=imdilate(kim2,strel('disk',1));
            end
            [cut2x,cut2y]=find(kim2 & perims);
            subcell = drawline(subcell,[cut1x(1) cut1y(1)],[cut2x(1) cut2y(1)],0);
     
            cutcell = cell;
            cutcell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4);    
            Lout(Lout==Lout(cutx,cuty))=0;
            Lout(cutcell==1) = chosencolor;
            for k = 2:max2(cutcell),
                Lout(cutcell==k) = max2(Lout)+k-1;
            end;
        else
            disp(['less than 2 perims! cell number: ' num2str(Lout(cutx,cuty)) ' in Lbot_back.'])
        end
    end
end
end
if OKorNot & ~finetuneimage
    Lout=renumberimage(Lout);
    for sfi=1:max2(Lout)
        [sfx,sfy]=find(Lout==sfi);
        if length(sfx)<100
            Lout(Lout==sfi)=0;
        end
    end
    Lout=renumberimage(Lout);
end
if pp(1) delete(pp(1));delete(pp(2));delete(pp(3));delete(pp(4));end;pp=0;
if bb delete(bb);delete(bbp);end;bb=0;
