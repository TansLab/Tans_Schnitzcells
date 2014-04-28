function theimg = MW_stampit(theimg,framenr)
%
%


global whitelist mywatermark; % MW

% EDIT by MW 7-3-2014
% This checks whether the frame is on the whitelist,
% if so, it adds a small circle ('watermark') to the plot
% as a flag to the user.   
%
if(inlist(whitelist,framenr))
    offset=5;
    for i = 1:size(mywatermark,1)
        for j = 1:size(mywatermark,2)
            if mywatermark(i,j,2)>0
                theimg(i+offset,j+offset,2) = 1;
            end
        end
    end
    %{
    figure(ourfig);
%         outpos=get(gcf,'Outerposition' );
    hold on
    h=imshow(mywatermark);
%         disp(num2str(size(mywatermark)));
%         set(gcf,'OuterPosition',outpos);
    set(h, 'AlphaData', mywatermark(:,:,2));
    hold off;
    %}
end

%
% END edit MW  