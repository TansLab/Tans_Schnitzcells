function SegImageConn= NW_ConnectDividedCellParts(SegImage,NegPhaseImage,AreaIdx)
%
% This function increases the area in 'Segimage' which has the number 'AreaIdx'
% until the Area touches a different segmented area. Prevents
% oversegmentation into subcell-areas.
% Function is called from NW_segphase_richMed.m
%
% In certain growth media (e.g. rich medium) cells have large internal
% contrast differences and one cell can be detected as several small units
% (typically 2-3. the poles are often seperated from the central part).
% This oversegmentation can be detected since the areas are suspiciously small
% for a normal cell.
% This function aims to close the gap between subareas by
% 1) Dilating the area of interes (AreaIdx)
% 2) Determining which of the newly added pixels (via dilation) have the
%    brightest values in a negative phase contrast (NegPhaseImage). 
%    Typically the missing pixels within cells in a negative phase image are
%    brighter than the connection to the neighbouring cells.
% 3) Adding the 20 (can be adjusted) brightest pixels to the area 'AreaIdx'
% 4) repeat the above steps until a connection to a neighbour cell is
%    established (overlap) or abort if a threshold number of rounds has been reached.
%
% *************************************
% INPUT
% *************************************
% SegImage:   Filled and labelled raw segmented image (before imclose correction), typically
%             a labelled C_seeds1  (bwlabel(..))
%             background must be =0
% NegPhaseImage: Inverted phase image, typically B_negPh
% AreaIdx:    Label of area which shall be increased
%
% ****************************************
% OUTPUT
%*****************************************
% SegImageConn: Cell 'AreaIdx' is increased to overlap a neighbouring cell.
%      The area label of this merged area is changed to the minimum number
%      of the numbers.



% ******** ADJUST **************
maxNumRounds=10;   % 5... max # of growth loops before abortion
Pxgrowth=20;         % 20... # of added pixels per loop

% ** structuring element **
% line along the long axis
cellangle=regionprops(SegImage==AreaIdx,'Orientation'); % major axis orientation.
cellangle=cellangle.Orientation;
mystruct_el=strel('line',4,cellangle); 
% alternative: disk: mystruct_el=strel('disk',2);   % disk,2.... structuring element for dilation
% *****************************

SegImageBackup=SegImage;

currentround=1;
overlapNeighbour=0;
SegImRest=(SegImage~=AreaIdx & SegImage~=0); % all cells, except for AreaIdx, are =1;

% ------------------------------------
% grow cell fragment
% ------------------------------------
while currentround<=maxNumRounds & overlapNeighbour==0
    
    SegImSub=(SegImage==AreaIdx); % only cell of interest =1. 
    if max2(SegImSub)==0
        disp('Cell has disappeared.')
        break
    end
   
    % dilate cell of interest
    SegImDil=imdilate(SegImSub,mystruct_el);    
    %newly added area after dilation    
    SegImDelta=SegImDil & ~ SegImSub; 
    % get pixel intensities of phase contrast image for delta area
    NegPhaseImageDelta=uint16(SegImDelta).* NegPhaseImage;
    
    % find the Pxgrowth (e.g. 20) largest pixel values and their location
    pxvalues=sort(NegPhaseImageDelta(:));
    pxvalues=pxvalues(pxvalues~=0);
    cutoffpxvalue=pxvalues(max((end-Pxgrowth),1)); % lowest accepted value
    %corresponding coordinates of the bright pixels:
    idxbrightpx=find(NegPhaseImageDelta>=cutoffpxvalue);
    
    %add new pixels to labelled seg-image (use the AreaIdx label, not '1')
    SegImage(idxbrightpx)=AreaIdx;
    
    currentround=currentround+1;
    % check overlap;
    %ovlerapNeighbour=sum(sum(SegImSub.*SegImRest));
    overlapNeighbour=sum(SegImRest(idxbrightpx));
    
end

% ------------------------------------
% include grown cell fragment into seg image if possible
% ------------------------------------
if overlapNeighbour==0 % max number of rounds exceeded
    SegImageConn=SegImageBackup;
    %disp(['Cell ' num2str(AreaIdx) ' could not be corrected because no neighbour found.'])
else % merge the newly touching cells
    SegImSub=(SegImage==AreaIdx); % only cell of interest =1.
    SegImCoverArea=SegImSub.*SegImageBackup; % only the area which is now covered by the newly grown cell is non-zero (has original labels)
    cellLabels=unique(SegImCoverArea);
    cellLabels=cellLabels(cellLabels~=0); % labels of covered/overlapping cell areas (ideally size=2);
    newLabel=min(cellLabels);
    % check if cell might be touching 2 other cells -> could lead to filled
    % background
    SegImSubDil=imdilate(SegImSub,strel('disk',1)); % enlargen by 1
    SegImDilCoverArea=SegImSubDil.*SegImageBackup; 
    possibletripleoverlap=unique(SegImDilCoverArea);
    possibletripleoverlap=possibletripleoverlap(possibletripleoverlap~=0);
    
    if length(cellLabels)~=2 | length(possibletripleoverlap)>=3;
       % disp(['Don''t know to which cell to connect to. Possible tripe overlap.' ...
       %     ' Will skip cell ' num2str(AreaIdx) '.'])
        SegImageConn=SegImageBackup;
    else
        for i=1:length(cellLabels)
            idx=cellLabels(i);
            SegImage(SegImage==idx)=newLabel;
        end
        SegImageConn=SegImage;
        %disp(['Corrected cell ' num2str(AreaIdx)]);
    end
end
    
    