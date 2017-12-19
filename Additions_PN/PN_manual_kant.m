function  [p,Lout,OKorNot,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,gotoframenum,DJK_settings] = ...
    PN_manual_kant(p,Lin,L_prec,phin,rect,rect_prec,phsub,DJK_settings,assistedCorrection)
% function  [p,Lout,OKorNot,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,gotoframenum,DJK_settings] = ...
%    PN_manual_kant(p,Lin,L_prec,phin,rect,rect_prec,phsub,DJK_settings,assistedCorrection)
%
% This function is called by PN_manualcheckseg.
%
% Input arguments:
% - p                   : holds general information about current movie
% - Lin,                : segmented input image
% - L_prec,             : segmented image of preceeding frame 
% - phin,               : phase image with this frame
% - rect,               : 
% - rect_prec,          : ROI recteangle of previous image
% - phsub,              : 
% - DJK_settings,       : 
% - assistedCorrection  : 
%
%
%
% output:
% crop_pop      obsolete and can be removed - MW

% *******
% TODO: include "updatedCellNumbers" in program (PN_imshowlabel) so that
% regionprops are only calculated for this numbers
% ********

% elapsed time for frame 444 in 2012-05-08. 390 cells

iptsetpref('imshowborder','tight');
backwards = 0;
global pos Limage ourfig pp phfig currentFrameStr figureToFocusOn showPhase

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


%% Main loop

while ~done
    clear j*  %j1=0;j2=0;
    
    %% Convert current framenumber to string
    if isfield(p,'currentFrame'), currentFrameStr = num2str(p.currentFrame);
    else currentFrameStr = '??'; warning('p.currentFrame not set'); end % MW TODO: num2str(p.currentFrame) can just be used directly below
     
    %% Update segmented figure
    
    % Newest version MW of showing figures             

    % Update the figure size
    if ~isfield(p, 'fullscreenResolution')
        
        % note imSize= leny, lenx, depth; screensize is lenx, leny
        imSize=size(phsub);
        desiredLocation = (p.screenSize([3,4])-[imSize([2,1])])/2;
        desiredPosition = [desiredLocation,imSize([2,1])];
        
        % left, bottom, width, height
        if ~any(desiredLocation+imSize([2,1])>p.screenSize([3,4]))
            set(ourfig, 'Position', desiredPosition); % to be given as x,y,width,height
            
            ha = axes('Parent',ourfig,'Units','pixels','Position',[1 1 imSize([2,1])]);
        else        
            % if figure is bigger than the screen, set it to the size of
            % the screen
            figure(ourfig); clf; % needs to be cleared to remove prev. img
            set(ourfig, 'position', p.screenSize,'Toolbar','none'); % hide toolbar because would cover img
            ha = axes('Parent',ourfig,'Units','pixels','Position',[1 1 p.screenSize(3:4)]);
        end        
    else
        % For full screen mode, just set both positions to full screen position
        figure(ourfig); clf; % needs to be cleared to remove prev. img
                                
        set(ourfig, 'units', 'normalized', 'outerposition', [0 0 1 1]); %p.fullscreenResolution+[FULLSCREENMARGIN FULLSCREENMARGIN -2*FULLSCREENMARGIN -2*FULLSCREENMARGIN],'Toolbar','none'); % hide toolbar because would cover img                        
        ha = axes('Parent',ourfig,'Units','normalized','Position',[0 0 1 1]);
        
        % old method
        %set(ourfig, 'position', p.fullscreenResolution,'Toolbar','none'); % hide toolbar because would cover img
        %ha = axes('Parent',ourfig,'Units','pixels','Position',[1 1 p.screenSize(3:4)]);
        
    end
        
    % Update the figure
    % NOTE: the figure needs to be cleared (matlab's clf function) to
    % prevent memory from filling up. In most cases this is done in the 
    % PN_imshowlabel function.
    figure(ourfig);
    if ~showPhase        
        if assistedCorrection && ~isempty(L_prec)        
            % Update figure with segmentation + suspicious cell detection
            PN_imshowlabel(p, Lout,rect,L_prec,rect_prec,'phaseImage',phsub);                        
        else
            % % Update figure with segmentation 
            PN_imshowlabel(p, Lout,[],[],[],'phaseImage',phsub);            
        end
    else
        % Show unsegmented image
        clf(gca); % to prevent memory filling up
        imshow(phsub,[]);
    end               
    
    set(ourfig,'WindowButtonMotionFcn', @MW_schnitzfigureinteraction);
    
    %iptsetpref('imshowborder','tight') % not sure if necessary
    
    % edited stuff
    %{
    figure(ourfig) % NW2012-05-10        
    
    % make the assisted image
    
    % If assisted correction desired call PN_imshowlabel with 
    % suspicious cell detection functionality
    if assistedCorrection && ~isempty(L_prec)       

        p.res=res;
        PN_imshowlabel(p, Lout,rect,L_prec,rect_prec,'phaseImage',phsub); %slow step 0.45sec! (NW 2012-05-10)
    else

        % If correction assistence not desired, call PN_imshowlabel
        % with [],[],[] values such that functionality is not invoked.
        % Note also here resized img is given as input, so resizing is
        % not necessary.
        if isfield(p,'res'), p=rmfield(p,'res'); end
        PN_imshowlabel(p,imresize_old(Lout,res),[],[],[],'phaseImage',imresize_old(phsub,res)); %slow step 0.2sec! (NW 2012-05-10)                    

    end                    
    
    % make the raw phase image
    if showPhase
        figure(phfig);
        imshow(phsub,[]);
        iptsetpref('imshowborder','tight'); % DJK 090111 added so Lc & phase overlap       
    end
    
    % Set function to interact with figure
    set(figureToFocusOn,'WindowButtonMotionFcn',@MW_schnitzfigureinteraction);      
    
    % Title of figure 
    set(ourfig,'name',['Frame ' currentFrameStr]);
        
    % set the sizes of those figures equal    
    if isfield(p, 'fullscreenResolution')
        % For full screen mode, just set both positions to full screen position
        set(ourfig, 'position', p.fullscreenResolution); 
        if showPhase, set(phfig, 'position', p.fullscreenResolution); end
    else
        % Normal mode (default)
        % calculate center position
        pos11                 = get(ourfig,'position'); % current position
        width = pos11(3); height = pos11(4); 
        x_left_bottom_screen  = 10 + (p.min_size(2)-width)/2 ;
        y_left_bottom_screen  = 35 + (p.min_size(1)-height)/2 ;

        % set figures to that position
        set(ourfig, 'position', [x_left_bottom_screen y_left_bottom_screen width height]); % DJK 090117          
        % also set phase image to that position
        if showPhase, set(phfig, 'position',[x_left_bottom_screen y_left_bottom_screen width height]); end

    end
    
    %theFigurePosition=get(ourfig,'Position');
    %set(phfig,'Position',theFigurePosition);
    
    % Set focus on desired figure (segmented img per default)    
    figure(figureToFocusOn); % MW TODO    
    %}
    
    % Old stuff
    %{
    set(0,'CurrentFigure',figureToFocusOn);    
    uistack(figureToFocusOn, 'top'); % required in matlab 2014 to actually put figure on top - MW
        % sometimes still does not work
        % --
        % * in that case figureToFocusOn is still set right
        % * pressing 'm' twice resolve the problem
    %}          
  
    
    %% Get keypress
    
    clear cc;
    validPress = 0;
    while ~validPress,
        
        ct=waitforbuttonpress; % wait for buttonpress
        
        % if it's a keyboard press get that character
        % (ct==1 when keyboard, 0 when mouse)
        if ct
            cc=get(gcf,'currentcharacter'); % obtain value button % MW 2015/01
            p.lastPressed = cc;
        else
            cc='*'; % flag for mouse click
        end
                        
        % Only count as valid if 
        if (ct==1) && isempty(cc), % if keyboard pressed (ct==1), check if valid key (not ctrl etc)
            validPress=0;
        else % ct==0: mouseclick
            validPress=1;
        end;                
        
    end;
    
    %% Analyze keyboard clicks
    
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
    if strcmp(cc,'m')
        
        % MW 2016/04
        % When upgrading to matlab 2014b some issues started appearing with
        % showing the images. Maybe displaying only the image of interest
        % is better?
        
        % MW 2015/01
        % Switching figure between phfig (phase image) and ourfig
        % (segmented image). Note that PN_manual_kant is terminated and 
        % called again when going to next img, so switch only remains valid
        % during editing of this frame.
        
        % And switch
        showPhase = ~showPhase;            
                
        if showPhase
            disp('Now showing phase only.');
        else
            disp('Now showing phase+segmentation.');
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
            if (p.showNr == 2) 
                if ~isfield(p,'slookup')
                    disp('No lookup table available to find schnitznrs (use MW_slookup.m to create p.slookup).');
                    p.showNr = 0;
                else
                    highestSchnitzIndx = max(p.slookup(size(p.slookup,1),:));
                    % Set up custom colormap
                    % easy way:
                    %theColorMap = linspecer(maxCellNo);
                    
                    % let's use one of the standard colormaps
                    standardColorMap = hsv(highestSchnitzIndx); % hsv jet
                    
                    % but mix it up such that neighbors have different
                    % colors                    
                    shuffle=randperm(highestSchnitzIndx); 
                    %shuffle = [10  20  28  35  29  27  64  40  33  37  47  58  22  36  18  21  50  57  34  25  19  43   1  49  16  60  23   3  48   9  45  38  44  46   4  26   7  15  54  59  55  24   6  14  42  13   8  61  63  41   5   2  30  53  31  52  32  51  56  17  11  62  12  39];
                    
                    % perform shuffling
                    standardColorMapShuffled = NaN(size(standardColorMap));
                    standardColorMapShuffled(:,1) = standardColorMap(shuffle,1);
                    standardColorMapShuffled(:,2) = standardColorMap(shuffle,2);
                    standardColorMapShuffled(:,3) = standardColorMap(shuffle,3);
                    
                    % how many copies of the colormap do we need?
                    %copiesNeeded = ceil(maxCellNo/COLORMAPSIZE);
                    
                    % create the color map
                    %theColorMap = repmat(standardColorMapShuffled,copiesNeeded,1);
                    
                    % set the colormap
                    p.customColors = [0 0 0; standardColorMapShuffled; 1 1 1];
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
            
            figure(addfig); % needed in case current figure changes (in Windows)
            imshow(makergb(+imresize_old(LZedge,5),imresize_old(Phzoom(:,:,1),5)));
            
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
            %subcolroi=imresize_old(~roipoly,1/res);
            subcolroi=~roipoly;
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
        elseif cc=='k' % cell joining option
            
            affectedcells = MW_joinmergers(p,Lout,rect,L_prec,rect_prec);            
            
            updatedCellNumbers=[updatedCellNumbers;affectedcells]; % note that this option is currently not used, ie. this code is useless
            
            % since L_prec has changed, let's reload it
            fileLocation=[p.segmentationDir p.movieName 'seg' sprintf('%03d', p.currentFrame-1) '.mat'];
            load(fileLocation,'Lc');
            L_prec=Lc;
            clear Lc; % just to prevent mixups
            
            done=0;
            
            %{
            disp(['obsolete and not correct version!'])
            %subcolroi=imresize_old(~roipoly,1/res);
            subcolroi=~roipoly;
            if size(subcolroi,1)~=size(Lout,1) | size(subcolroi,2)~=size(Lout,2)
                subcolroi2=zeros(size(Lout));
                subcolroi2(1:size(subcolroi,1),1:size(subcolroi,2))=subcolroi;
                subcolroi=~subcolroi2;
            end
            Lout_undo=Lout;   %for undo step NW 2014-01
            Lout=(double(Lout).*double(subcolroi));
            %           Lout=renumberimage(Lout);
            
            done=0;
            %}
            
        elseif cc=='t'
            
            %subcolroi=imresize_old(~roipoly,1/res);
            subcolroi=~roipoly;
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
            
            %subcolroi=imresize_old(~roipoly,1/res);
            subcolroi=~roipoly;
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
        %elseif cc == 's' % MW option was double used?
        %    savetemp=1;
        %    done=1;
        elseif cc == 'n' 
            
            % Addition MW to go fullscreen
            % Note that you need to go to next frame after toggling this
            % option off.
            
            % obtain monitor resolution
            % NOTE: this is not needed when using normalized units, which
            % is done now, so this can be removed at some point.. (TODO)
            set(0,'units','pixels');            
            Pix_SS = get(0,'screensize');
            
            % Set flag to go fullscreen
            if ~isfield(p,'fullscreenResolution')
                p.fullscreenResolution = Pix_SS;
            else
                p = rmfield(p,'fullscreenResolution');
            end           
            
            % figure needs to be cleared because otherwise smaller
            % version of figure will remain visible
            figure(ourfig); clf;
            
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
        elseif cc == 'r' % 
            if ~isfield(p,'showmubar')
                p.showmubar=1;
                disp('Show mu bar activated');
            else
                p.showmubar=~p.showmubar;
                disp('Show mu bar toggled.');
            end                                   
            
            % Old code for previous function
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
                  
                [Lout, Lout_undo] = MW_removestuffnotoverlappingwithprevious(p,Lout,L_prec,rect,rect_prec);
                
            else
                warning('Previous segmentation image not loaded or not existent.')
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
                cutcell_dilated = carefuldilate( cutcell, strel('diamond',1), 2, Lcell); 
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
                    set(ourfig,'WindowButtonMotionFcn', @MW_schnitzfigureinteraction);
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
                

                % perim is perimeter of dilated cell
                perim = bwperim(imdilate(subcell,strel('disk',1)));
                

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
        
        % left mouse button: connect two cells
        % ===
        if cz(1)=='n'
            Lout_undo=Lout;   %for undo step NW 2014-01
            
            % join cells
            pos1=pos;
            j1=Lout(round(pos(1,2)),round(pos(1,1)));
            figure(ourfig);
            set(ourfig,'WindowButtonMotionFcn',@MW_schnitzfigureinteraction);
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
        
            done=0;
            
        % .............
        % ===
        elseif (cz(1)=='e' & Lout(round(pos(1,2)),round(pos(1,1))))
            Lout_undo=Lout;   %for undo step NW 2014-01
            % erase cell
            Lout(Lout==Lout(round(pos(1,2)),round(pos(1,1))))=0;
            
            % alternate) = ctrl+left button or right mouse button
            % Will cut at this location
            
            % updated cellnumbers NW2012-05-10
            updatedCellNumbers=[updatedCellNumbers;Lout(round(pos(1,2)),round(pos(1,1)))];
            
            done=0;
            
        % Right mouse click: cut cell where clicked
        % ===
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
            
            
            % perim is perimeter of dilated cell
            perim = bwperim(imdilate(subcell,strel('disk',1)));
            
            
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
            
            done=0;
            
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