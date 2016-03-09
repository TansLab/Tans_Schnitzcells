%by Philippe Nghe 16/01/3012
function cutPoints = PN_FindBranchPoints(skelim)

cutPoints = zeros(size(skelim));

%find points from which 3 branches emerge
dskelim = double(skelim);
ibranch = imfilter(dskelim,ones(3));
ibranch(~skelim) = 0;
brchIm = ismember(ibranch,4);

%go through the branching zones
cc=bwconncomp(brchIm);
stats = regionprops(cc,'Centroid');
for ii = 1:cc.NumObjects
    
    %create a local zoom tool
    centerP = int16(stats(ii).Centroid);
    imDisk = zeros(size(skelim));
    imDisk(centerP(2),centerP(1))=1;
    imDisk = imdilate(imDisk,strel('disk',6));  
    
    %if the 3 branches are long enough, only cut the "orthogonal" one, else savagely remove the branching zone
    localSkelet = imDisk & skelim;
    ccLocal = bwconncomp(imDisk & ~skelim,4);
    if ccLocal.NumObjects == 3 
        statLocal = regionprops(ccLocal,'Area');
        [~,idx] = max([statLocal.Area]); %detect the biggest "solid" angle
        mainFragment = ismember(labelmatrix(ccLocal), idx);
        mainFragment = bwmorph(mainFragment,'dilate'); %ceate a protection area of the skelet
        cutPoints(localSkelet & ~mainFragment) = 1;
    else
        cutPoints(localSkelet) = 1;
    end

end

clear dskelim ibranch brchIm cc stats centerP imDisk localSkelet ccLocal statLocal mainFragment
end




