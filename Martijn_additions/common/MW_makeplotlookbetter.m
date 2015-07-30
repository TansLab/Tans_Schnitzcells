
function MW_makeplotlookbetter(FontSize)
% function MW_makeplotlookbetter(FontSize)
%
% Currently only sets all font sizes to <FontSize>. 
% Maybe in future add some more conveniences.

%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',FontSize,'fontWeight','normal')
set(gca,'FontSize',FontSize)

end