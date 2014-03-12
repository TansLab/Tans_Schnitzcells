% To plot all schnitz raw data together as listed in the file AllData.
% Plots should be the same as the ones used in GetMultSwTimes2.

figure
hold all
clear legendstring
for i=1:size(AllData,1)
    
%     styleMap=['r-';'b-';'g-';'k-';'m-';'c-'];
%     style=char(styleMap(1+mod(i,size(styleMap,1))));
      colors=jet(size(AllData,1))
      morecolors=spring(size(AllData,1))
      
      Tbr=AllData{i,4}(1);
      Tsw=AllData{i,5}(1);
 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     FOR FYmins vs MYs :

%     %to plot based on their real time, first parent, then the switching schnitz :
%     plot([AllData{i,2}{1,1} AllData{i,2}{1,end}],[AllData{i,3}{1,1} AllData{i,3}{1,end}],style);
    

%     % to make the switching & parental plots lying on the same starting point:
%     plot([AllData{i,2}{1,1}-AllData{i,2}{1,1}(1) AllData{i,2}{1,end}-AllData{i,2}{1,1}(1)],...
%          [AllData{i,3}{1,1} AllData{i,3}{1,end}],style);
    
     % to synchronise the division previous to switching lies at the same time point:
       % for late switching schnitz, the daughter schnitz is included and the
       % end of the switching schnitz is no longer clear.
%      if Tsw>Tbr
%      plot([AllData{i,2}{1,1}-AllData{i,2}{1,end}(1) AllData{i,2}{1,end}-AllData{i,2}{1,end}(1)],...
%           [AllData{i,3}{1,1} AllData{i,3}{1,end}],'Color',colors(i,:));
%      else plot([AllData{i,2}{1,1}-AllData{i,2}{1,1}(1) AllData{i,2}{1,end}-AllData{i,2}{1,1}(1)],...
%           [AllData{i,3}{1,1} AllData{i,3}{1,end}],'Color',morecolors(i,:));
%      end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FOR phase vs MYs :  

%      %only for switching schnitz:
%      plot([AllData{i,11}{1,end}],...
%      [AllData{i,3}{1,end}]);
     
%      %for both parent & switching schnitz, but discontinuous:
%      plot([AllData{i,11}{1,1}-1 AllData{i,11}{1,end}],...
%      [AllData{i,3}{1,1} AllData{i,3}{1,end}]);
 

% %    TO PLOT ALSO THE REGRESSION fot MYs vs synchronised time: 
%      x1=AllData{i,2}{1,1}
%         y1=(AllData{i,13}(1)*x1+AllData{i,13}(2))
%      plot(x1-AllData{i,2}{1,end}(1),y1,'-r');
%     
%      x2=AllData{i,2}{1,end}
%         y2=(AllData{i,14}(1)*x2+AllData{i,14}(2))
%      plot(x2-AllData{i,2}{1,end}(1),y2,'-r');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
legendstring(i)=AllData{i,1}(1);
end
legend(legendstring,'Location','Best')

hold off


