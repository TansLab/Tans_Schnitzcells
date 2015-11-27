function MW_showtwoframeslabeled(LcFrame1, LcFrame2,h)

    % PARAMETERS
    FONTSIZE=10; 
    FIGUREHEIGHT=1000;
    FIGURESCREENMARGIN=40;
    halfFontSize = round(FONTSIZE/2);
        
    % Get region properties frame n and n+1 (resp. 1 and 2)
    propL1 = regionprops(LcFrame1,'Centroid'); % area
    propL2 = regionprops(LcFrame2,'Centroid'); % area

    % Get size of frame 1
    [sizex, sizey] = size(LcFrame2)        
    
    % Was necessary when using axis (old)
    % set(gca,'YDir','Reverse')

    % pad frame 1
    padLcFrame1 = padarray(LcFrame1,size(LcFrame2)-size(LcFrame1),'post');

    % modify matrices for display
    frame1Pic = (1-(padLcFrame1>0).*.3);
    frame2Pic = (1-(LcFrame2>0).*.5);
    imshow(frame2Pic.*frame1Pic)

    % Labels for frame 1 (n)
    for i = 1:numel(propL1)

        % position of label
        textx=propL1(i).Centroid(1)-halfFontSize;
        texty=propL1(i).Centroid(2)-halfFontSize;        
        
        % print label
        text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[0,0,1],'FontWeight','bold')            
    end


    % Labels for frame 2 (n+1)
    for i = 1:numel(propL2)

        textx=propL2(i).Centroid(1)-halfFontSize;
        texty=propL2(i).Centroid(2)-halfFontSize;

        text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[0,0,0],'FontWeight','bold')
    end

    % Set figure size
    %set(h,'Position',[1,1,sizex,sizey])
    set(h,'Position',[FIGURESCREENMARGIN,FIGURESCREENMARGIN,FIGUREHEIGHT,ceil(sizey/sizex*FIGUREHEIGHT)])
    
end