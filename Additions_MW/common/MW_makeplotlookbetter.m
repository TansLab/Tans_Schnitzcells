
function MW_makeplotlookbetter(FontSize,optionalParameters,plotSize,setRenderer)
% function MW_makeplotlookbetter(FontSize,optionalParameters,plotSize,setRenderer)
%
% Currently only sets all font sizes to <FontSize>. 
%
% Input arguments:
% - FontSize        
%                   Determines font size that is set
% - optionalParameters.style
%                   Will cause some additional changes to be made based on
%                   this string.
% - plotSize        if set will set plotsize to 2*plotSize=2*[width height]; 
%                   note that the size is doubled to increase resolution
%                   for non-vector images; you should also double the
%                   font-size manually accordingly in the input arguments.
% - setRenderer     if setRenderer=1 will set renderer to painters
%
% One can use function "linspecer" to generate colors that are
% distinguishable by human eye.
%
% Example:
% MW_makeplotlookbetter(8*2,[],SIZE,1)

if ~exist('optionalParameters','var')
    optionalParameters.style = 'default';
end

%Set all fontsizes
set(findall(gcf,'type','text'),'FontSize',FontSize,'fontWeight','normal','FontName','Arial')
set(gca,'FontSize',FontSize)
set(gca,'FontName','Arial')

if exist('plotSize','var')
    set(gcf,'Units','centimeters','Position',[3,3,plotSize*2]);
    
    % Now for pdfs the following is also important:
    %set(gcf, 'PaperUnits', 'normalized')
    %set(gcf, 'PaperPosition', [0 0 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', plotSize*2)
    set(gcf, 'PaperPosition', [0 0 plotSize*2])
end

if exist('setRenderer','var')
    set(gcf,'RendererMode','manual','Renderer','Painters');
end



%{
if strcmp(optionalParameters.style, 'CBmanuscript')
    FIGSIZE=[7.5 5.63]; OFFSET=[5 5];
    set(hLenLifetime, 'units', 'centimeters', 'pos', [OFFSET OFFSET+FIGSIZE*2]);
end
%}

end




