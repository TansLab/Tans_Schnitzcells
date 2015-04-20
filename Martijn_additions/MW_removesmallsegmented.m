function L = MW_removesmallsegmented(L, minimum)
% MW 2015/04
% subfunction that removes small areas (# pixels < minimum) in segmented
% image matrix L. L contains indices, which tells you to which segmented
% cell the pixel / matrix entry belongs to.
% 
% bwareaopen could also be used, but cannot distinguish cells. (And as such
% cannot handle areas that are bordering tightly.)

% Loop over all known cell indices
for index = [1:max(L(:))]   
    % Find all pixels thought to belong to that cell
    hits = L==index;
    % Calculate total amount of pixels
    pixelcount = sum(hits(:));
    % Remove if there's too little of them
    if pixelcount < minimum
        L(find(hits)) = 0;
    end
end

end