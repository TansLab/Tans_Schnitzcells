function theimg = MW_stampit(theimg,framenr)
%
%
% Function to show user which frames have already been approved. It will
% add a green circle to frames that are in list of approved frames. This
% list is called "whitelist".
%
% TODO MW: this function is somewhat redundant, because whitelist is not
% necessary to determine whether a frame is approved. (See e.g. also
% DJK_analyzeseg.m.)

% For efficiency, the watermark image is already loaded in
% MW_manualcheckseg.m.

global whitelist mywatermark; % MW

% This checks whether the frame is on the whitelist,
% if so, it adds a small circle ('watermark') to the plot
% as a flag to the user.   
if(ismember(framenr,whitelist))
    offset=5;
    for i = 1:size(mywatermark,1)
        for j = 1:size(mywatermark,2)
            if mywatermark(i,j,2)>0
                theimg(i+offset,j+offset,2) = 1;
            end
        end
    end
    
    % Old way of doing it:
    %{
    figure(ourfig);
    hold on
    h=imshow(mywatermark);
    set(h, 'AlphaData', mywatermark(:,:,2));
    hold off;
    %}
    
end