
function [generations] =  distancelca(p, schnitzcells, schnitzIdx1, schnitzIdx2)
    % function [generations] =  givemeschnitzinfoforframe(p, schnitzcells, schnitzIdx1, schnitzIdx2)
    %
    % This functions returns the distance between the last common ancestor
    % (lca), defined as the sum of the number of generations from two given
    % schnitzes to the lca, divided by two.
    
    % Settings
    MAXITERATIONS = 100;
    
    totalDistance = NaN;
    
    index1 = schnitzIdx1;
    index2 = schnitzIdx2;
    
    ancestry1 = [];
    ancestry2 = [];
    
    lcaFound = 0; deadend = 0; iterations = 0;
    while ~lcaFound && ~deadend && ~(iterations>MAXITERATIONS)
       
        % expand the list of ancestry on generation up
        ancestry1(end+1) = index1;
        ancestry2(end+1) = index2;
        
        % look if the list share values of entries
        [hits1, distancestolca1] = ismember(ancestry1, ancestry2);
        [hits2, distancestolca2] = ismember(ancestry2, ancestry1);
        
        % if there are hits the lca is found
        if ~isempty(find(hits1)) && ~isempty(find(hits2)) % TODO, could be programmed better?
            totalDistance = (find(distancetolca1)-1+find(distancetolca2)-1)/2; % -1 bc nr of linkes, not elements
                % note that I assume there is only one distance found; 
                % (i.e. one non-zero value in distancestolca2). this should be the case
                
            lcaFound = 1; % output found
        else
            % is there are parent for schnitz 1? look at that one
            candidateParent1 = schnitzcells(index1).P;
            if ~((isempty(candidateParent1) || candidateParent1==0))                
                index1 = candidateParent1;
            end
            % idem for schnitz 2
            candidateParent2 = schnitzcells(index2).P;
            if ~((isempty(candidateParent2) || candidateParent2==0))                
                index2 = candidateParent2;
            end
            
            % if there's no ancestry left for both; give up; no lca found
            % (output will be NaN)
            if ((isempty(candidateParent1) || candidateParent1==0)) && ((isempty(candidateParent2) || candidateParent2==0))
                deadend = true;
            end
        end
        
        iterations=iterations+1;
        
    end

    generations = totalDistance;
    
    if iterations>MAXITERATIONS
        disp('OOPS! Maximum iterations reached... This is not good.');
    end

end