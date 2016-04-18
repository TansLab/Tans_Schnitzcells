function linklistschnitz = MW_convertsimplelinklisttoschnitzlist(linklist, debugmode);
% function linklistschnitz = MW_convertsimplelinklisttoschnitzlist(linklist, debugmode);
%
% Converts simple list linking cellno's between two frames to schnitzcell's
% format for doing this.
% INPUT
% - linklist:  simple link list 
%           [cellno_fr1, cellno_fr2; cellno_fr1, cellno_fr2; etc.]
% - debugmode: 
%           set to 0. Set to 1 to debug.

%%
linklistschnitz = [linklist(:,1), zeros(size(linklist,1),2), linklist(:,2)];
uniqueWithoutZeros=unique(linklistschnitz(:,1));
uniqueWithoutZeros=uniqueWithoutZeros(find(uniqueWithoutZeros>0));
[parentcount,indices] = hist(linklistschnitz(:,1),uniqueWithoutZeros);

parentlist = indices(find(parentcount>1))';
cellListFrame1 = linklist(:,1);
dividesInLinkList = find(changem_clone(cellListFrame1,parentcount,indices)>1);

if debugmode
    dividingrowsbefore = linklistschnitz(dividesInLinkList,:)
end

movecount = ones(1,max(indices))+1;
for idx=dividesInLinkList'
    oldline = linklistschnitz(idx,:);
    newline = [0, 0, 0, linklistschnitz(idx,4)];
    parentNumber = linklistschnitz(idx,1);
    if movecount(parentNumber)<4
        newline(movecount(parentNumber))=parentNumber;
        movecount(parentNumber)=movecount(parentNumber)+1;
        linklistschnitz(idx,:) = newline;
    else
        disp('ERROR: parent has more than 2 daughters! Leaving >2 untouched!');
    end
end

% sort
linklistschnitz = sortrows(linklistschnitz,4);
    
if debugmode
    dividingrowsafter = linklistschnitz(dividesInLinkList,:)    
end

end