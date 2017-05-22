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
% expand 2xN linklist to 4xN
linklistschnitz = [linklist(:,1), zeros(size(linklist,1),2), linklist(:,2)];
% unique values in left column 
uniqueWithoutZeros=unique(linklistschnitz(:,1));
uniqueWithoutZeros=uniqueWithoutZeros(uniqueWithoutZeros>0);
% now count the nr of instances of each unique entry
indices = uniqueWithoutZeros';
binedges = [indices-.5 indices(end)+.5];
parentcount = histcounts(linklistschnitz(:,1),binedges);

% identify which ones are occuring multiple times (cells that have more
% than one cells linked have probably divided).
parentlist = indices(find(parentcount>1))';
% get all cells linked in frame 1
cellListFrame1 = linklist(:,1);

% identify the cells that have divided
% (create a list where instead of cell indices their count is now listed)
dividesInLinkList = find(changem_clone(cellListFrame1,parentcount,indices)>1);
if debugmode
    dividingrowsbefore = linklistschnitz(dividesInLinkList,:)
end

% go over cells that divide
movecount = ones(1,max(indices))+1; % which column to place parent 
for idx=dividesInLinkList'
    % retrieve old link entry
    oldline = linklistschnitz(idx,:); 
    % start making new one
    newline = [0, 0, 0, linklistschnitz(idx,4)]; 
    % get parent nr
    parentNumber = linklistschnitz(idx,1); 
    % while spot 2 and 3 are available
    if movecount(parentNumber)<4 
        % fill 'm up with the parent nr
        newline(movecount(parentNumber))=parentNumber; 
        % then make sure next entry this parent other spot is taken
        movecount(parentNumber)=movecount(parentNumber)+1;
        % put in new line
        linklistschnitz(idx,:) = newline;
    else
        % If a parent has >2 children, make the >2 ones barren.
        newline = [0, 0, 0, linklistschnitz(idx,4)]; 
        linklistschnitz(idx,:) = newline;
        % But give error
        warning('ERROR: parent has more than 2 daughters! Making >2 barren!');
    end
end

% sort
linklistschnitz = sortrows(linklistschnitz,4);
    
if debugmode
    dividingrowsafter = linklistschnitz(dividesInLinkList,:)    
end

end