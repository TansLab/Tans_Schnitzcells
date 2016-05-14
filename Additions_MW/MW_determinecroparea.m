function [leftTop, rightBottom] = MW_determinecroparea(p, range)
% 'little helper' function determinecroparea(p, range)
% ===
% Tries to determine last image in range, and simply calls getrect for that
% to then output leftTop and rightBottom positions for cropping.

    % Determine and load last phase image
    lastIndex = range(end)+(range(end)-range(end-1)); % otherwise issue going from prelim. to full analysis
    fileNameSearchString = [p.imageDir, [p.movieName,'-p-*' num2str(lastIndex) '.tif'] ]
    DlastPhase = dir(fileNameSearchString);
    if numel(DlastPhase)==0
        disp(['Looking for ' fileNameSearchString]);
        warning(['No images found w. index ' num2str(lastIndex)]);

        % Same code as above, but for other lastIndex
        warning(['Since margin failed, just using last img from range.']);
        lastIndex = range(end); 
        fileNameSearchString = [p.imageDir, [p.movieName,'-p-*' num2str(lastIndex) '.tif'] ];
        DlastPhase = dir(fileNameSearchString);
        
        if numel(DlastPhase)==0
            disp(['Looking for ' fileNameSearchString]);
            disp(['Also no images found w. index ' num2str(lastIndex)]);
            error('Also this failed.. Wrong path?');
        end
    end
    imToLoad = DlastPhase(round(numel(DlastPhase)/2));
    myImg = imread([p.imageDir imToLoad.name]);
    myImgSize = size(myImg);
    
    if (range(end)-range(end-1))>1
        disp('I''m taking some extra frames after the last in range since your deltaframe>1.');
    end

    % Show figure and getrect
    h=figure(); imshow(myImg,[]);
    myCrop = round(getrect());

    % Get leftTop and rightBottom of cropping region
    leftTop = myCrop(1:2)
    rightBottom = myCrop(1:2)+ myCrop(3:4)

    % Correct for even/uneven requirements
    if mod(leftTop(1),2) == 0 % if even
        leftTop(1) =  leftTop(1) -1 ; % make uneven
    end
    if mod(leftTop(2),2) == 0 % if even
        leftTop(2) =  leftTop(2) -1 ; % make uneven
    end
    if mod(rightBottom(1),2) == 1 % if uneven
        rightBottom(1) =  rightBottom(1) +1 ; % make even
        if rightBottom(1) > myImgSize(2), rightBottom(1) = rightBottom(1) - 2; end % make sure not to expand beyond img size
    end
    if mod(rightBottom(2),2) == 1 % if uneven
        rightBottom(2) =  rightBottom(2) +1 ; % make even
        if rightBottom(2) > myImgSize(1), rightBottom(2) = rightBottom(2) - 2; end % make sure not to expand beyond img size
    end        
    % output to user    
    disp(['leftTop=' mat2str(leftTop) '; rightBottom=' mat2str(rightBottom)])

    % close figure
    close(h);

end