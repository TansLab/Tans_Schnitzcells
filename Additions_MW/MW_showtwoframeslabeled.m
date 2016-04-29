function MW_showtwoframeslabeled(LcFrame1, LcFrame2,h,showlabels)
%function MW_showtwoframeslabeled(LcFrame1, LcFrame2,h,showlabels)

    %% PARAMETERS
    FONTSIZE=10; 
    %FIGUREHEIGHT=1000;
    %FIGURESCREENMARGIN=40;
    halfFontSize = round(FONTSIZE/2);        
    if ~exist('showlabels','var'), showlabels=1;  end
    
    % Get region properties frame n and n+1 (resp. 1 and 2)
    propL1 = regionprops(LcFrame1,'Centroid'); % area
    propL2 = regionprops(LcFrame2,'Centroid'); % area

    % Get size of frame 1
    [sizex, sizey] = size(LcFrame2);
    
    % Was necessary when using axis (old)
    % set(gca,'YDir','Reverse')

    % Padding such that frames are equal size
    howMuchPadding=[0, 0; 0, 0]; % [padding fr1i, padding fr1j, padding fr2i, padding fr2j]
    whichOneBiggerBothDirections = size(LcFrame2)<=size(LcFrame1);
    whichToChange = whichOneBiggerBothDirections+1;
    howMuchPadding(whichToChange(1),1) = abs(size(LcFrame2,1)-size(LcFrame1,1));
    howMuchPadding(whichToChange(2),2) = abs(size(LcFrame2,2)-size(LcFrame1,2));
    padLcFrame1 = padarray(LcFrame1,howMuchPadding(1,:),'post');
    padLcFrame2 = padarray(LcFrame2,howMuchPadding(2,:),'post');    
        
    % modify matrices for display
    frame1Pic = (1-(padLcFrame1>0).*.3);
    frame2Pic = (1-(padLcFrame2>0).*.5);
    % colored version
    frame1Pic = (1-(padLcFrame1>0)*.5);
    frame2Pic = (1-(padLcFrame2>0)*.5);
    frame1PicColored = ones([size(frame1Pic),3]);
    frame1PicColored(:,:,1) = frame1Pic; % color blue channel
    frame1PicColored(:,:,2) = frame1Pic; % color blue channel
    frame2PicColored = ones([size(frame2Pic),3]);
    frame2PicColored(:,:,2) = frame2Pic; % color red channel
    frame2PicColored(:,:,3) = frame2Pic; % color red channel
    % show them
    %imshow(frame1Pic.*frame2Pic)
    imshow(frame1PicColored.*frame2PicColored)

    %% Labels for frame 1 (n)
    if showlabels
        for i = 1:numel(propL1)

            % position of label
            textx=propL1(i).Centroid(1)-halfFontSize;
            texty=propL1(i).Centroid(2)-halfFontSize;        

            % print label
            text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[0,0,1],'FontWeight','bold');            
        end
    end

    %% Labels for frame 2 (n+1)
    if showlabels
        for i = 1:numel(propL2)

            textx=propL2(i).Centroid(1)-halfFontSize;
            texty=propL2(i).Centroid(2)-halfFontSize;

            text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[1,0,0],'FontWeight','bold');
        end
    end

    % Set figure size
    %set(h,'Position',[1,1,sizex,sizey])
    %set(h,'Position',[FIGURESCREENMARGIN,FIGURESCREENMARGIN,FIGUREHEIGHT,ceil(sizey/sizex*FIGUREHEIGHT)]);
    set(h,'units','normalized','outerposition',[0 0 1 1])
    
end