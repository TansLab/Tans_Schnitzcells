function theimg = MW_stampit(theimg,p)
% function theimg = MW_stampit(theimg,p)
%
%
% Function to show user which frames have already been approved. It will
% add a green circle to frames that are in list of approved frames.
%
% Whether the current frame is approved is saved in 
% p.CurrentFrameApprovedFlag.
%
% For efficiency, the watermark image is already loaded in
% MW_manualcheckseg.m. It is stored in p.mywatermark.
%
% function theimg = MW_stampit(theimg,p)

% If a frame is approved, a small circle ('watermark') is added to the plot
% as a flag to the user.   
if isfield(p, 'CurrentFrameApprovedFlag')
if p.CurrentFrameApprovedFlag    
    
    offset=5;
    for i = 1:size(p.mywatermark,1)
        for j = 1:size(p.mywatermark,2)
            if p.mywatermark(i,j,2)>0 % if green channel watermark value >0
                theimg(i+offset,j+offset,2) = 1; % set green channel img
            end
        end
    end
end

%figure, imshow(theimg,[]) TODO REMOVE

end