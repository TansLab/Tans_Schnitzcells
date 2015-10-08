function  [p,Lout,OKorNot,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,gotoframenum,DJK_settings] = ...
    PN_manual_kant(p,Lin,L_prec,phin,rect,rect_prec,phsub,DJK_settings,assistedCorrection)
% function  [p,Lout,OKorNot,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,gotoframenum,DJK_settings] = ...
%    PN_manual_kant(p,Lin,L_prec,phin,rect,rect_prec,phsub,DJK_settings,assistedCorrection)
%
% Input arguments:
% - p                   : holds general information about current movie
% - Lin,                : segmented input image
% - L_prec,             : segmented image of preceeding frame 
% - phin,               : phase image with this frame
% - rect,               : 
% - rect_prec,          : 
% - phsub,              : 
% - DJK_settings,       : 
% - assistedCorrection  : 
%
%
%

% *******
% TODO: include "updatedCellNumbers" in program (PN_imshowlabel) so that
% regionprops are only calculated for this numbers
% ********

% elapsed time for frame 444 in 2012-05-08. 390 cells



iptsetpref('imshowborder','tight');
backwards = 0;
global pos Limage ourfig res pp phfig currentFrameStr

%set(phfig,'Visible', 'on') % ugly way to force phase image away from foreground during computation time NW2012-10-05


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

% create array with cell numbers that changed (-> regioprops have to be
% recalculated) (NW2012-05-10)
updatedCellNumbers=[];
%

% initialize image for undo-step ('u')
Lout_undo=Lout;

figureToFocusOn = ourfig;

while ~done
    clear j*  %j1=0;j2=0;
    
    % Convert current framenumber to string
    if isfield(p,'currentFrame'), currentFrameStr = num2str(p.currentFrame);
    else currentFrameStr = '??'; warning('p.currentFrame not set'); end % MW TODO: num2str(p.currentFrame) can just be used directly below
    
    % ** NW 2012-05-10 **
    
    if ~ishandle(ourfig) % if segmented image figure doesn't exist yet (should never be the case)
       
        disp('WARNING: unexpected behavior detected in PN_manual_kant, segmented image figure doesn''t exist.');
        
        figure(ourfig);
        clf reset;
    
        if assistedCorrection && ~isempty(L_prec)
            Lshow = PN_imshowlabel(p,Lout,rect,L_prec,rect_prec,'phaseImage',phsub);
            Lshow = imresize_old(Lshow,res);
            imshow(Lshow);
            
        else
            PN_imshowlabel(p,imresize_old(Lout,res),[],[],[],'phaseImage',imresize_old(phsub,res)); % MW, [],[],[] signals no susp. cell. detect.
        end
    
        pos11 = get(phfig,'position'); % current position
        set(ourfig, 'position', pos11); % DJK 090117
    
         set(ourfig,'name',['Pos: ',num2str(pos(1,2)),' , ',num2str(pos(1,1)),...
            '  Val: ',num2str(double(Limage(pos(1,2),pos(1,1))))]);
    else        
        
        % Update segmented figure
        figure(ourfig) % NW2012-05-10        
        
        if assistedCorrection && ~isempty(L_prec)
            
            % If assisted correction desired call PN_imshowlabel with 
            % suspicious cell detection functionality
            
            p.res=res;
            PN_imshowlabel(p, Lout,rect,L_prec,rect_prec,'phaseImage',phsub); %slow step 0.45sec! (NW 2012-05-10)
            
            %{
            % MW REMOVED 2015/07
            
            Lshow = PN_imshowlabel(p, Lout,rect,L_prec,rect_prec,'phaseImage',phsub); %slow step 0.45sec! (NW 2012-05-10)
            %PN_imshow.. is slow step (0.45 sec)
            
            % resize and show
            Lshow = imresize_old(Lshow,res);                        
            imshow(Lshow,'InitialMagnification','fit');
            %}
            
        else
            
            % If correction assistence not desired, call PN_imshowlabel
            % with 0,0,0 values such that functionality is not invoked.
            
            % clf reset; %NW2012-05-10 <> NOT CLEAR DESIRED??!
            PN_imshowlabel(p,imresize_old(Lout,res),[],[],[],'phaseImage',imresize_old(phsub,res)); %slow step 0.2sec! (NW 2012-05-10)            
            % Note that for historic reasons PN_imshowlabel calls imshow
            % when 3-5th parameters are [],[],[].
            % Note also here resized img is given as input, so resizing is
            % not necessary.
            
            % MW 2015/01
            % set(ourfig, 'position', pos11); % DJK 090117 % a little awkward to redo, but should only affect first image
            %stop4b=toc
        end                    
                
        %  sometimes setting the image name leads to an error
        % -> create dummy position value
        % NW2014-04
        if pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2)
             curr_val1=num2str(double(Limage(pos(1,2),pos(1,1))));
        else
             curr_val1=-1;
        end
        
        % Set phase and segmented image to same position
        pos11 = get(phfig,'position');         
        set(ourfig, 'position', pos11); 
        
        % Set to fullscreen if desired
        if isfield(p, 'fullscreenResolution')
            set(ourfig, 'position', p.fullscreenResolution); 
            set(phfig, 'position', p.fullscreenResolution);             
        end
        
        % Title of figure 
        set(ourfig,'name',['Frame ' currentFrameStr ', Pos: ',num2str(pos(1,2)),' , ',num2str(pos(1,1)),...
                '  Val: ',num2str(curr_val1)]);
            
        % Set focus on desired figure (segmented img per default)
        set(0,'CurrentFigure',figureToFocusOn);

    end                    
    
    % Define function to retrieve mouseclick position.
    set(ourfig,'WindowButtonMotionFcn',['global pos currentFrameStr Limage ourfig res pp phfig;pos=max(1,round((1/res)*get(gca,''CurrentPoint'')));',...
        'if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2));',...
        'curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));else;curr_val=''-1'';end;',...
        'set(ourfig,''name'',[''Frame '' currentFrameStr '', Pos: '',num2str(pos(1,2)),'' , '',num2str(pos(1,1)),',...
        '''  Val: '',curr_val]);']);
    %   '''  Val: '',curr_val]);pp=showbox(phfig,pos,pp,size(Limage),res);figure(ourfig);']);
    
    % MW 2015/01 also enable for phase fig (strictly not necessary)
    set(phfig,'WindowButtonMotionFcn',['global pos currentFrameStr Limage ourfig res pp phfig;pos=max(1,round((1/res)*get(gca,''CurrentPoint'')));',...
        'if (pos(1,2)>0 & pos(1,2)<size(Limage,1) & pos(1,1)>0 & pos(1,1)<size(Limage,2));',...
        'curr_val=num2str(double(Limage(pos(1,2),pos(1,1))));else;curr_val=''-1'';end;',...
        'set(phfig,''name'',[''Frame '' currentFrameStr '', Pos: '',num2str(pos(1,2)),'' , '',num2str(pos(1,1)),',...
        '''  Val: '',curr_val]);']);
        
    % blubb
    % here the regionprops of the next and previous image could already be
    % calculated to quicken update of image when frame changes
    % NW 2012-05-10
        
    validPress = 0;
    while ~validPress,
        
        ct=waitforbuttonpress; % wait for buttonpress
        cc=get(gcf,'currentcharacter'); % obtain value button % MW 2015/01
        
        if (ct==1) && isempty(cc), % if keyboard pressed (ct==1), check if valid key (not ctrl etc)
            validPress=0;
        else % ct==0: mouseclick
            validPress=1;
        end;
        
    end;
    
    % set(ourfig,'WindowButtonMotionFcn',''); MW REMOVE THIS LINE? Don't
    % see point since above Fn was set. MW 2015/01
    
    %if ct cc=get(1,'currentcharacter');else cc=get(1,'selectiontype');end
    
    % ****************************************************
    % *********** START KEYBOARD CLICKS ******************
    % ****************************************************
    
    % bring phase image to front (and background again)
    % treated independently of 'real' actions since update in segImage
    % colors is not wanted. [NW 2014-04]
    
    
 % 'm' switch between phase contrast and segmentation image
    if ct & cc=='m'
        
        % MW 2015/01
        % Switching figure between phfig (phase image) and ourfig
        % (segmented image). Note that PN_manual_kant is terminated and 
        % called again when going to next img, so switch only remains valid
        % during editing of this frame.
        
        % And switch
        if figureToFocusOn == phfig % if current is ph 
            figureToFocusOn = ourfig; % switch to segmented
        elseif figureToFocusOn == ourfig % vice versa
            figureToFocusOn = phfig;
        else % if other figure is on focus (not to be the case)
            figureToFocusOn = ourfig; % switch to default, i.e. segmented one
        end
            
        done = 0; % continue keypress-loop
        
    end           
                
    if ct 

        % cc=get(ourfig,'currentcharacter');
           
        % save, continue etc
        if cc==' '
            OKorNot=1;
            done=1;
        elseif cc=='q'
            Lout=Lin;
            OKorNot=0;
            quit_now=1;
            done=1;
        elseif cc=='s' % MW 2015/07 addition
            % enable/disable (schnitzcells) numbering
            % showNr options:
            % 0) Don't show numbers
            % 1) Show label numbers (numbering from segmentation)
            % 2) Show schnitz numbers            
            if ~isfield(p,'showNr')
                % showNr wasn't active yet, will activate to show nrs in general
                p.showNr = 1;
            else
                % change the options
                p.showNr = p.showNr+1; % 
            end
            if p.showNr > 2, p.showNr = 0; end
            if (p.showNr == 2) && ~isfield(p,'slookup')
                disp('No lookup table available to find schnitznrs (use MW_slookup.m to create p.slookup).');
                p.showNr = 0;
            end            
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
     % *** deactivated (anyways never worked) (NW2015-07, 'b' is used for a different action now) 
     %   elseif cc=='b' %does not work.
     %       if bb delete(bb);delete(bbp);end;bb=0;
     %       if bbs
     %           bbs=0;
     %       else
     %           cutx=round(pos(1,2));
     %           cuty=round(pos(1,1));
     %           if Lout(cutx,cuty)
     %              [fx,fy] = find(Lout==Lout(cutx,cuty));
     %               figure(phfig);
     %               bb=plot(fy*res,fx*res,'w.');
     %               set(bb,'markersize',2);
     %               bbp=plot(cuty*res,cutx*res,'r.');
     %               set(bbp,'markersize',24);
     %               figure(ourfig);
     %               bbs=1;
     %           end
     %       end
     % *** 
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
            imshow(makergb(+imresize_old(LZedge,5),imresize_old(Phzoom(:,:,1),5)));
            figure(addfig); % needed in case current figure changes (in Windows)
            subaddcell=imresize_old(roipoly,1/5);%(phin);
            if max2(Lzoom(subaddcell>0))>0
                disp('overlaps existing cell; ignored.');
            else
                Lout_undo=Lout;   %for undo step NW 2014-01
                Lzoom(subaddcell)=max2(Lout)+1;
                Lout(zoomrect(1):zoomrect(2),zoomrect(3):zoomrect(4))=Lzoom;
                updatedCellNumbers=[updatedCellNumbers;max2(Lout)+1]; % blubb
            end
            close(addfig)
            figure(ourfig)
            done=0;
        elseif cc=='x'
            Lout_undo=Lout;   %for undo step NW 2014-01
            subcolroi=imresize_old(~roipoly,1/res);
            if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
                subcolroi2=zeros(size(Lout));
                subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
                subcolroi=subcolroi2;
            end
            
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout=(double(Lout).*double(subcolroi));
            %           Lout=renumberimage(Lout); %don't use if you want to
            %           work with "updatedCellNumbers"!
            % get cellnumbers that were affected by blacking out
            affectedcells=unique(Lout.*(1-subcolroi));
            affectedcells=affectedcells(affectedcells>0);
            updatedCellNumbers=[updatedCellNumbers;affectedcells]; %blubb
            %
            done=0;
        elseif cc=='k' % don't use!
            disp(['obsolete and not correct version!'])
            subcolroi=imresize_old(~roipoly,1/res);
            if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
                subcolroi2=zeros(size(Lout));
                subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
                subcolroi=~subcolroi2;
            end
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout=(double(Lout).*double(subcolroi));
            %           Lout=renumberimage(Lout);
            
            done=0;
            
        elseif cc=='t'
            
            subcolroi=imresize_old(~roipoly,1/res);
            if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
                subcolroi2=zeros(size(Lout));
                subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
                subcolroi=subcolroi2;
            end
            [xe,ye]=find(edge(subcolroi)==1);
            Lout_undo=Lout;   %for undo step NW 2014-01
            for qe=1:length(xe)
                if Lout(xe(qe),ye(qe))~=0
                    Lout(Lout==Lout(xe(qe),ye(qe)))=0;
                end
            end
            Lout=(double(Lout).*double(subcolroi));
            %               Lout=renumberimage(Lout);
            % get cellnumbers that were affected by blacking out
            affectedcells=unique(Lout.*(1-subcolroi));
            affectedcells=affectedcells(affectedcells>0);
            updatedCellNumbers=[updatedCellNumbers;affectedcells]; %blubb
            %
            done=0;
            
        elseif cc=='v' %NW2014-02 opposite of terace 't' -> restricts to region of interest
            
            subcolroi=imresize_old(~roipoly,1/res);
            if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
                subcolroi2=zeros(size(Lout));
                subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
                subcolroi=subcolroi2;
            end
            [xe,ye]=find(edge(subcolroi)==1);
            Lout_undo=Lout;   %for undo step NW 2014-01
            for qe=1:length(xe) % examine edge of roi
                if Lout(xe(qe),ye(qe))~=0
                    Lout(Lout==Lout(xe(qe),ye(qe)))=0;
                end
            end
            Lout=(double(Lout).*double(~subcolroi)); %inverted to 't'. % examine area of roi
            %               Lout=renumberimage(Lout);
            % get cellnumbers that were affected by blacking out
            affectedcells=unique(Lout.*subcolroi); %TOBE CHECKED! unique(Lout.*(1-subcolroi));
            affectedcells=affectedcells(affectedcells>0);
            updatedCellNumbers=[updatedCellNumbers;affectedcells]; %blubb
            %
            done=0;
            
        %elseif cc == 'c'
        %    disp(['don''t use!'])
        %    crop_pop=1;
        %    done=1;
        %    extra= 15; % <- extra number of pixels on either side of highlighted region
        %    LNfull=zeros(p.fullsize(1),p.fullsize(2));
        %    LNfull(rect(1):rect(3), rect(2):rect(4))=Lout;
        %    [fx,fy]= find(LNfull);
        %    xmin= max(min(fx) - extra, 1);
        %    xmax= min(max(fx) + extra, size(LNfull,1));
        %    ymin= max(min(fy) - extra, 1);
        %    ymax= min(max(fy) + extra, size(LNfull,2));
        %    newrect= [xmin ymin xmax ymax];
        elseif cc == 's'
            savetemp=1;
            done=1;
        elseif cc == 'n' 
            
            % Addition MW to go fullscreen
            % Note that you need to go to next frame after toggling this
            % option off.
            
            % obtain monitor resolution
            set(0,'units','pixels');            
            Pix_SS = get(0,'screensize');
            % Set flag to go fullscreen
            if ~isfield(p,'fullscreenResolution')
                p.fullscreenResolution = Pix_SS;
            else
                p = rmfield(p,'fullscreenResolution');
            end
            
            % (Setting the size here doesn't work since it's edited later.)
            %set(ourfig, 'units','pixels','position',Pix_SS);  % handles are also stored in phfig or ourfig
            %set(ourfig, 'units','normalized','position',[0 0 1 1]);  % handles are also stored in phfig or ourfig
        elseif cc == 'e'
            disp('"expand image" currently not in use.');
            %         expand_segimage(p,imnum,expandvalue);
            %         OKorNot=1;
            %         done=1;
            %         dontsave=1;
            %         backwards = 2;
        elseif cc == 'f'
            DJK_settings.finetuneimage = 1;
        elseif cc == 'g'
            gotoframenum = input('goto loopindex = ');
            OKorNot=1;
            done=1;
            dontsave=1;
            % blubb! special case because everything has to be
            % recalculated!!! NW 2012-05-10. updatedCellNumbers=ALL
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
            disp(['currently not in use']) %NW 2012-05-10
            %   Lout=renumberimage(Lout); 
        elseif cc == 'r' % 
            disp(['currently not in use']) %NW 2012-05-10
             %   % renumber this cell
             %   cutx=round(pos(1,2));
             %   cuty=round(pos(1,1));
             %   chosencolor=Lout(cutx,cuty);
             %   cell=zeros(size(Lout));
             %   cell(Lout==Lout(cutx,cuty))=Lout(Lout==chosencolor);
             %   [fx,fy] = find(cell);
             %   xmin = max(min(fx)-5,1);
             %   xmax = min(max(fx)+5,size(cell,1));
             %   ymin = max(min(fy)-5,1);
             %   ymax = min(max(fy)+5,size(cell,2));
             %   subcell = cell(xmin:xmax, ymin:ymax);
             %   cell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4);
             %   Lout(Lout==Lout(cutx,cuty))=0;
             %   Lout(cell==1) = chosencolor;
             %   for k = 2:max2(cell),
             %       Lout(cell==k) = max2(Lout)+k-1;
             %   end;
        elseif cc == 'l' % MW 2014/12
            % Switch to show perimiter of bacteria instead of area
            if ~existfield(p, 'showPerim')
                % if not there create switch
                p.showPerim = 1; % and turn it on
            else
                p.showPerim = ~p.showPerim; % otherwise just switch it
            end
        elseif cc == 'o'
            % obliterate all but this cell
            cutx=round(pos(1,2));
            cuty=round(pos(1,1));
            chosencolor=Lout(cutx,cuty);
            % basically all cells have been updated (=deleted)
            updatedCellNumbers=unique(Lout);
            %
            if (chosencolor > 0)
                Lout_undo=Lout;   %for undo step NW 2014-01
                Lout = (Lout==chosencolor);
            end
        elseif cc==27
            Lout=Lin;OKorNot=0;
            done=1;
        elseif cc =='',
            disp('you typed shift');
        elseif cc=='i' & Lout(round(pos(1,2)),round(pos(1,1)))         % Noreen 2012-05. fills 1 cell       
            cellonlyImage=(Lout==Lout(pos(1,2),pos(1,1))); % =1 for chosen cell, =0 elsewhere
            cellonlyImage=imfill(cellonlyImage,'holes');
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout(cellonlyImage==1)=Lout(pos(1,2),pos(1,1));            
            % very crude: update all cells (NW2012-10-05)
            updatedCellNumbers=unique(Lout);
            %
        elseif cc=='h'                              % Philippe 2012-02
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout = PN_reseed(Lout,phsub,round(pos(1,2)),round(pos(1,1)));
            % BLUBB very crude and probably only the highest cell number or
            % merged cell number has to be updated. NW2012-05-10
            updatedCellNumbers=unique(Lout);
            %
        elseif cc=='j'
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout=imfill(Lout,'holes');
        elseif cc=='c'    %NW 2014-01 morphologicaly close each cell area.
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout=NW_imclose_eachArea(Lout,strel('disk',10));
        elseif cc=='d' %           remove possible dirt. Dirt often has a large (NW2012-10)
            Lout_undo=Lout;   %for undo step NW 2014-01
             for runcell=1:max2(Lout)      % convex hull compared to its actual area.
                 % update (NW 2013-01): dirt has large perimeter(outline)
                 % compared to area. More reliable than convex hull (bent
                 % cells!)
                 % update (NW 2013-01-30): combine properties of convex
                 % hull and perimeter. further excludes all areas with
                 % holes (fill ells first!)
                 deletecellConv=0;
                 deletecellPerim=0;
                 deletecellHole=0;
                 areaprops=regionprops(Lout==runcell,'Area','ConvexArea','Perimeter','EulerNumber');
                 allarea=[areaprops.Area];
                 allconvex=[areaprops.ConvexArea];
                 allperim=[areaprops.Perimeter];
                 containsholes=[areaprops.EulerNumber]; % (holes: value<=0)
                 
                 ratioconv=sum(allarea)/sum(allconvex);
                 if ratioconv<0.82%0.88 %blubb . maybe adjust (around 0.85. dirt=low values)
                     %                      large value= excludes a lot of areas (maybe cells))
                     deletecellConv=1;
                 end                 
                 ratioperim=sum(allperim)/sum(allarea);
                 if ratioperim>0.17%0.17 %blubb . maybe adjust (around 0.17. dirt=high values)
                               %                      small value= excludes a lot of areas (maybe cells))
                     deletecellPerim=1;
                 end
                 if containsholes<=0
                     deletecellHole=1;
                 end
                 
                 if deletecellHole==1  | (deletecellConv==1 &  deletecellPerim==1) % maybe change to 'or'/'and' condition
                     Lout(Lout==runcell)=0;
                 end
             end
        elseif cc=='u'   % undo one step (NW 2014-01)
            Lout=Lout_undo;
                     
       elseif cc=='b'   % restrict segmentation to overlaps with previous segmentation image (useful for non-growing cells)
                        % NW2015-07
            if ~isempty(L_prec) % previous image exists (can only be loaded in asisted correction mode)
                Lout_undo=Lout;   %for undo step NW 2014-01
                % old and new seg need to be located properly:
                Lout_full=zeros(p.fullsize);
                Lout_full(rect(1):rect(3),rect(2):rect(4))=Lout;
                L_prec_full=zeros(p.fullsize);
                L_prec_full(rect_prec(1):rect_prec(3),rect_prec(2):rect_prec(4))=L_prec;
                
                L_prec_mask=L_prec_full;             
                %L_prec_mask=imdilate(L_prec_full,strel('disk',3)); % enlargen prev. image overlap mask if needed
                L_prec_mask=L_prec_mask>0;  % binary mask from previous image

                cellsall=unique(Lout);  % same numbers as Lout_full
                cellsoverlap=unique(Lout_full.*L_prec_mask); %cells overlapping with previous seg. image
                cellstodelete=setdiff(cellsall,cellsoverlap);
                for i=cellstodelete(:)'
                    Lout(Lout==i)=0;
                end
            else
                disp('Previous segmentation image not loaded or not existent.')
            end
            done=0; 
                        
                        
            %%%%%%%%%%%%%%%%%%% START DJK 090105 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif cc=='9'  % still part of keyboard clicks
            if DJK_settings.fill_cut
                DJK_settings.fill_cut = 0;
                disp('Filling of Cutting Pixels Turned OFF');
            else
                DJK_settings.fill_cut = 1;
                disp('Filling of Cutting Pixels Turned ON');
            end
        elseif cc=='0'
            if DJK_settings.figs
                DJK_settings.figs = 0;
                disp('Showing Figures Turned OFF');
            else
                DJK_settings.figs = 1;
                disp('Showing Figures Turned ON');
            end
            
        elseif ( (cc=='1' || cc=='2' || cc=='3') && Lout(round(pos(1,2)),round(pos(1,1))) )
            % 1: New Phase Cut
            % 2: More original, but better hard cut with random dilate afterwards
            % 3: Original cut of line between closest perim
            
            Lout_undo=Lout;   %for undo step NW 2014-01
            
            % cutx & cuty are coordinate that is clicked
            cutx = round(pos(1,2)); cuty = round(pos(1,1)); % disp(['clicked on (' num2str(cutx) ',' num2str(cuty) ')']);
            
            % color of clicked cell
            chosencolor = Lout(cutx,cuty);
            
            % get only the cell that is gonna be cut
            Lcell = (Lout == chosencolor);
            
            % cut the cell in the selected fashion, returns subcell with correct colors
            if cc=='1'
                cutcell = DJK_manual_cutPhase(Lcell, phsub, cutx, cuty, DJK_settings.figs);
            elseif cc=='2'
                %still working on it: cutcell = DJK_manual_cut(Lcell, cutx, cuty, DJK_settings.figs);
            elseif cc=='3'
                cutcell = DJK_manual_cutOriginal(Lcell, cutx, cuty, DJK_settings.figs);
            end
            
            % remove cells that are just small pixels
            r = regionprops(cutcell, 'area');
            fpts = find([r.Area]<5);
            for i = 1:length(fpts)
                cutcell(cutcell == fpts(i))= 0;
            end
            cutcell = DJK_randomRenumberImage(cutcell);
            disp(['rarely used function! probably update of renumbered cells is wrong!!!'])
            
            % dilate into original cell
            if (DJK_settings.fill_cut)
                cutcell_dilated = carefuldilate( cutcell, strel('diamond',1), 2, Lcell); %figure(111);imshowlabel(cutcell);figure(112);imshowlabel(cutcell_dilated);pause;close(111);close(112);
                cutcell = cutcell_dilated; % not required with carefuldilate: cutcell(find(Lcell>0)) = cutcell_dilated(find(Lcell>0));
            end
            
            % remove original cell
            Lout(Lcell) = 0;
            
            % place cutcell (can be 1 or more cells)
            cellnos = unique(cutcell); % cellnos(1) is background
            Lout(find(cutcell == cellnos(2)))= chosencolor; % first cell has original color
            for j = 3:length(cellnos)
                Lout(find(cutcell == cellnos(j)))= max2(Lout)+1;
            end
            
            disp(['# cells now : ' str3(length(unique(Lout))-1)]);
            
            % rarely used functions. for security, impose to recalc
            % regionproperties of all cells NW2012-05-10
            updatedCellNumbers=unique(Lout);
            %%%%%%%%%%%%%%%%%%%  END DJK 090105  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % ------------ start NW 2013-05 replace mouse ---------
        elseif (cc=='4' || cc=='7' || cc=='5')
            % alternative ways to replace mouse button function (does the same things as mouse button functions)
                if cc=='4' % join cells by pressing '4' twice
                    Lout_undo=Lout;   %for undo step NW 2014-01
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
                    cc=get(ourfig,'currentcharacter');
                    set(ourfig,'WindowButtonMotionFcn','');
                    if ct==0 & cc=='4'
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

                    % updated cellnumbers NW2012-05-10
                    updatedCellNumbers=[updatedCellNumbers;j1;j2];

            % remove cell = shift+left button (or both?) -> press '7'to remove
            % cell
            elseif (cc=='7' & Lout(round(pos(1,2)),round(pos(1,1))))
                
                Lout_undo=Lout;   %for undo step NW 2014-01
                % erase cell
                Lout(Lout==Lout(round(pos(1,2)),round(pos(1,1))))=0;

                % alternate) = ctrl+left button or right mouse button
                % Will cut at this location

                % updated cellnumbers NW2012-05-10
                updatedCellNumbers=[updatedCellNumbers;Lout(round(pos(1,2)),round(pos(1,1)))];

            elseif (cc=='5' & Lout(round(pos(1,2)),round(pos(1,1)))) %-> press '5' to cut cell into 2
                Lout_undo=Lout;   %for undo step NW 2014-01

                % cutx & cuty are coordinate that is clicked
                cutx = round(pos(1,2));
                cuty = round(pos(1,1));

                % extract cell that will be cut
                chosencolor = Lout(cutx,cuty);
                cell = zeros(size(Lout));
                cell(Lout==Lout(cutx,cuty)) = Lout(Lout==chosencolor);
                [fx,fy] = find(cell);
                xmin = max(min(fx)-5,1);
                xmax = min(max(fx)+5,size(cell,1));
                ymin = max(min(fy)-5,1);
                ymax = min(max(fy)+5,size(cell,2));
                subcell = cell(xmin:xmax, ymin:ymax);
                % subcell is only cell that will be cut
                % figure(1);imshowlabel(subcell);pause;close(1);

                % perim is perimeter of dilated cell
                perim = bwperim(imdilate(subcell,strel('disk',1)));
                % figure(1);imshowlabel(perim);pause;close(1);

                % starting from clicked point, will increase a box untill 2 sides are
                % found: this will be perims
                perims = zeros(size(perim));
                radp = 1;
                while max2(perims)<2 & radp<41
                    pxmin = max(cutx-xmin+1-radp,1);
                    pxmax = min(cutx-xmin+1+radp,size(perims,1));
                    pymin = max(cuty-ymin+1-radp,1);
                    pymax = min(cuty-ymin+1+radp,size(perims,2));
                    perims(pxmin:pxmax,pymin:pymax) = bwlabel(perim(pxmin:pxmax,pymin:pymax));
                    radp = radp+1;
                end
                % figure(1);imshowlabel(perims);pause;close(1);

                % if indeed 2 sides are found, will cut
                if max2(perims)>1
                    % kim is image with only clicked point drawn
                    kim=zeros(size(subcell));
                    kim(cutx-xmin+1,cuty-ymin+1)=1;

                    % look for start of drawline
                    kim1=kim;
                    % increase size of kim untill it hits perims
                    while ~any(any(kim1 & perims))
                        kim1=imdilate(kim1,strel('disk',1));
                    end
                    % randomly select first point as start of drawline
                    [cut1x,cut1y]=find(kim1 & perims);

                    % now go for end of drawline, first remove points of side of start from perims
                    color1=perims(cut1x(1),cut1y(1));
                    perims(perims==color1)=0;
                    kim2=kim;
                    while ~any(any(kim2 & perims))
                        kim2=imdilate(kim2,strel('disk',1));
                    end
                    % randomly select first point as end of drawline
                    [cut2x,cut2y]=find(kim2 & perims);
                    color2=perims(cut2x(1),cut2y(1));

                    % cut cell by drawing line
                    subcell = drawline(subcell,[cut1x(1) cut1y(1)],[cut2x(1) cut2y(1)],0);

                    % cutcell is original cell, but now with seperate cells different colors
                    cutcell = cell;
                    cutcell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4);

                    % remove original cell in Lout
                    Lout(Lout==Lout(cutx,cuty))=0;

                    % first cell gets original color
                    Lout(cutcell==1) = chosencolor;



                    % new cells get new color
                    for k = 2:max2(cutcell),

                        % updated cellnumbers NW2012-05-10
                        updatedCellNumbers=[updatedCellNumbers;chosencolor; max2(Lout)+k-1];

                        Lout(cutcell==k) = max2(Lout)+k-1;
                    end;



                else
                    disp(['less than 2 perims! cell number: ' num2str(Lout(cutx,cuty)) ' in Lbot_back.'])
                end
                end % end replace mouse clicks (if cc='4')

            % ----------- end NW 2013-05 replace mouse button functions --------
            
        end
    
    % ****************************************************
    % *********** START MOUSE KLICKS *********************
    % ****************************************************    
        
    else 

        cz=get(ourfig,'selectiontype');
        
        % normal = left mouse button
        if cz(1)=='n'
            Lout_undo=Lout;   %for undo step NW 2014-01
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
            
            % updated cellnumbers NW2012-05-10
            updatedCellNumbers=[updatedCellNumbers;j1;j2];
            
            % delete cell = shift+left button (or both?)
        elseif (cz(1)=='e' & Lout(round(pos(1,2)),round(pos(1,1))))
            Lout_undo=Lout;   %for undo step NW 2014-01
            % erase cell
            Lout(Lout==Lout(round(pos(1,2)),round(pos(1,1))))=0;
            
            % alternate) = ctrl+left button or right mouse button
            % Will cut at this location
            
            % updated cellnumbers NW2012-05-10
            updatedCellNumbers=[updatedCellNumbers;Lout(round(pos(1,2)),round(pos(1,1)))];
            
        elseif (cz(1)=='a' & Lout(round(pos(1,2)),round(pos(1,1))))
            Lout_undo=Lout;   %for undo step NW 2014-01
                        
            % cutx & cuty are coordinate that is clicked
            cutx = round(pos(1,2));
            cuty = round(pos(1,1));
            
            % extract cell that will be cut
            chosencolor = Lout(cutx,cuty);
            cell = zeros(size(Lout));
            cell(Lout==Lout(cutx,cuty)) = Lout(Lout==chosencolor);
            [fx,fy] = find(cell);
            xmin = max(min(fx)-5,1);
            xmax = min(max(fx)+5,size(cell,1));
            ymin = max(min(fy)-5,1);
            ymax = min(max(fy)+5,size(cell,2));
            subcell = cell(xmin:xmax, ymin:ymax);
            % subcell is only cell that will be cut
            % figure(1);imshowlabel(subcell);pause;close(1);
            
            % perim is perimeter of dilated cell
            perim = bwperim(imdilate(subcell,strel('disk',1)));
            % figure(1);imshowlabel(perim);pause;close(1);
            
            % starting from clicked point, will increase a box untill 2 sides are
            % found: this will be perims
            perims = zeros(size(perim));
            radp = 1;
            while max2(perims)<2 & radp<41
                pxmin = max(cutx-xmin+1-radp,1);
                pxmax = min(cutx-xmin+1+radp,size(perims,1));
                pymin = max(cuty-ymin+1-radp,1);
                pymax = min(cuty-ymin+1+radp,size(perims,2));
                perims(pxmin:pxmax,pymin:pymax) = bwlabel(perim(pxmin:pxmax,pymin:pymax));
                radp = radp+1;
            end
            % figure(1);imshowlabel(perims);pause;close(1);
            
            % if indeed 2 sides are found, will cut
            if max2(perims)>1
                % kim is image with only clicked point drawn
                kim=zeros(size(subcell));
                kim(cutx-xmin+1,cuty-ymin+1)=1;
                
                % look for start of drawline
                kim1=kim;
                % increase size of kim untill it hits perims
                while ~any(any(kim1 & perims))
                    kim1=imdilate(kim1,strel('disk',1));
                end
                % randomly select first point as start of drawline
                [cut1x,cut1y]=find(kim1 & perims);
                
                % now go for end of drawline, first remove points of side of start from perims
                color1=perims(cut1x(1),cut1y(1));
                perims(perims==color1)=0;
                kim2=kim;
                while ~any(any(kim2 & perims))
                    kim2=imdilate(kim2,strel('disk',1));
                end
                % randomly select first point as end of drawline
                [cut2x,cut2y]=find(kim2 & perims);
                color2=perims(cut2x(1),cut2y(1));
                
                % cut cell by drawing line
                subcell = drawline(subcell,[cut1x(1) cut1y(1)],[cut2x(1) cut2y(1)],0);
                
                % cutcell is original cell, but now with seperate cells different colors
                cutcell = cell;
                cutcell(xmin:xmax,ymin:ymax) = bwlabel(subcell,4);
                
                % remove original cell in Lout
                Lout(Lout==Lout(cutx,cuty))=0;
                
                % first cell gets original color
                Lout(cutcell==1) = chosencolor;
                
                
                
                % new cells get new color
                for k = 2:max2(cutcell),
                    
                    % updated cellnumbers NW2012-05-10
                    updatedCellNumbers=[updatedCellNumbers;chosencolor; max2(Lout)+k-1];
                    
                    Lout(cutcell==k) = max2(Lout)+k-1;
                end;
                
                
                
            else
                disp(['less than 2 perims! cell number: ' num2str(Lout(cutx,cuty)) ' in Lbot_back.'])
            end
        end
    end
    
    % get correct array for updated cell numbers NW2012-05-10. Not yet used
    updatedCellNumbers=unique(updatedCellNumbers);
    updatedCellNumbers=updatedCellNumbers(updatedCellNumbers>0);
    %
end
if ~DJK_settings.finetuneimage & OKorNot 
    % TODO!! if you want to use updatedNumberCells
     %tic %start time
    Lout=NW_renumberimage(Lout);
     %toc  % elapsed time: 0.03 sec
    
     %for sfi=1:max2(Lout) %slow
     %   [sfx,sfy]=find(Lout==sfi);
     %   if length(sfx)<100
     %       Lout(Lout==sfi)=0;
     %   end
     %end
     propArea=regionprops(Lout,'Area'); %faster
     miniCells = find([propArea.Area] < 100);
     for i=miniCells
         Lout(Lout==i)=0;
     end
     
    Lout=NW_renumberimage(Lout);
    %toc % elapsed time: 0.07sec
end
if ~DJK_settings.finetuneimage & (savetemp==2 & done==1) %case 'w' (NW2012-05-10)
    % renumber but don't delete small cell fragments
    % TODO!! if you want to use updatedNumberCells
    Lout=NW_renumberimage(Lout);
end


if pp(1) delete(pp(1));delete(pp(2));delete(pp(3));delete(pp(4));end;pp=0;
if bb delete(bb);delete(bbp);end;bb=0;