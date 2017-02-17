
function MW_makeplotlookbetter(FontSize,optionalParameters)
% function MW_makeplotlookbetter(FontSize)
%
% Currently only sets all font sizes to <FontSize>. 
%
% Input arguments:
% - FontSize        
%                   Determines font size that is set
% - optionalParameters.style
%                   Will cause some additional changes to be made based on
%                   this string.
%
% One can use function "linspecer" to generate colors that are
% distinguishable by human eye.

if ~exist('optionalParameters','var')
    optionalParameters.style = 'default';
end

%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',FontSize,'fontWeight','normal','FontName','Arial')
set(gca,'FontSize',FontSize)
set(gca,'FontName','Arial')

%{
if strcmp(optionalParameters.style, 'CBmanuscript')
    FIGSIZE=[7.5 5.63]; OFFSET=[5 5];
    set(hLenLifetime, 'units', 'centimeters', 'pos', [OFFSET OFFSET+FIGSIZE*2]);
end
%}

end




