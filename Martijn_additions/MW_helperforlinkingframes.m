

FRAMEN=155

FONTSIZE=12; 
halfFontSize = round(FONTSIZE/2);

filename1 = [p.segmentationDir,p.movieName,'seg',str3(FRAMEN),'.mat'];
load(filename1);
LcFrame1 = Lc;
filename2 = [p.segmentationDir,p.movieName,'seg',str3(FRAMEN+1),'.mat'];
load(filename2);
LcFrame2 = Lc;


propL1 = regionprops(LcFrame1,'Centroid'); % area
propL2 = regionprops(LcFrame2,'Centroid'); % area

[sizex, sizey] = size(LcFrame2)

h = figure(1); clf; hold on;
set(gca,'YDir','Reverse')
%plot(0);
xlim([0,sizex]);
ylim([0,sizey]);

% Labels for frame N
for i = 1:numel(propL1)
    l1=plot(propL1(i).Centroid(1),propL1(i).Centroid(2),'.','Color',[0,0,1])
    
    textx=propL1(i).Centroid(1)-halfFontSize;
    texty=propL1(i).Centroid(2)-halfFontSize;
    %if assistedCorrection, textx=textx*p.res; texty=texty*p.res; end    
    text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[0,0,1],'FontWeight','bold')            
end


% Labels for frame N+1
for i = 1:numel(propL2)
    l2=plot(propL2(i).Centroid(1),propL2(i).Centroid(2),'.','Color',[1,0,0])
    
    textx=propL2(i).Centroid(1)-halfFontSize;
    texty=propL2(i).Centroid(2)-halfFontSize;
    %if assistedCorrection, textx=textx*p.res; texty=texty*p.res; end
    
    text(textx,texty,sprintf('%03d', i),'FontSize',FONTSIZE,'Color',[1,0,0],'FontWeight','bold')
end

legend([l1,l2],{'N','N+1'})

xlabel('pixels');
ylabel('pixels');
