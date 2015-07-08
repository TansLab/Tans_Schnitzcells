function linklistschnitz = MW_linkframes(p, frameNumber)
% function linklistschnitz = MW_linkframes(p, frameNumber)
%
% Performs tracking for frames frameNumber and frameNumber+1.
%
% if p.debug is valid field, then also figures are plotted.

%% Parameters 

DISKSIZE=15;
MARGIN=10;

if ~exist('frameNumber')
    frameNumber=190;
end

if isfield(p,'debug')
    debugmode=1;
else
    debugmode=0;
end

%% Load data
data=load(['F:\A_Tans1_step1_incoming_not_backed_up\2015-06-02\pos7crop\segmentation\pos7cropseg'  sprintf('%03d', frameNumber) '.mat'],'LNsub');
frame1=data.LNsub;
if debugmode
    figure(1)
    PN_imshowlabel(p,frame1,[],[],[]);
end

data=load(['F:\A_Tans1_step1_incoming_not_backed_up\2015-06-02\pos7crop\segmentation\pos7cropseg' sprintf('%03d', frameNumber+1) '.mat'],'LNsub');
frame2=data.LNsub;
if debugmode
    figure(2), PN_imshowlabel(p,frame2,[],[],[]);
end

myColors = [0,0,0; distinguishable_colors(max([frame1(:); frame2(:)]),[0,0,0])];

%%
% DISKSIZE defined at top

% Resize such that they are equally big
size(frame1)
size(frame2)

se = strel('disk',DISKSIZE,8);

% dilate frame 1
frame1BW = (frame1>0);
frame1BWarea=imdilate(frame1BW,se);
frame1BWarea=imerode(frame1BWarea,se);
% dilate frame 2
frame2BW = (frame2>0);
frame2BWarea=imdilate(frame2BW,se);
frame2BWarea=imerode(frame2BWarea,se);

%figure(3), imshow(frame2BWwork)


%% Perimiter image
% Find 
perimImg = bwperim(frame2BWarea);
%perimImg = imdilate(perimImg,strel('disk',3,4)); % imdilate use TODO can be optimized
areaImg = frame2BWarea+frame2;

if debugmode
    figure(4)
    PN_imshowlabel(p,areaImg,[],[],[]);
end

%%
STATS1 = regionprops(frame1BWarea, 'Area', 'Centroid', 'EulerNumber');
STATS2 = regionprops(frame2BWarea, 'Area', 'BoundingBox', 'Centroid', 'EulerNumber');
    
if (numel(STATS1) > 1) || (numel(STATS2) > 1)
    error('Failed to detect single colony area.');
end

ratio = STATS2.Area/STATS1.Area
    
frame1resized = imresize(frame1, ratio);
if debugmode
    figure(6), imshow(frame1resized)
end
   
%% create aligned frame #1
% MARGIN defined at top

% Create frame that has equal size to frame 2, but contains info from
% frame1, with centroid at same position as frame 2.

frame1recentered = zeros(size(frame2));
deltacentroids = round(STATS2.Centroid-STATS1.Centroid);
sizeframe1 = size(frame1);
frame1recentered([1:sizeframe1(1)]+deltacentroids(1), ...
                 [1:sizeframe1(2)]+deltacentroids(2)) = frame1;

             
             
%% plot overlay of aligned images
             
if debugmode
    overlaycentered = frame1recentered+frame2;
    figure, PN_imshowlabel(p,overlaycentered,[],[],[]);
end

%% multiply the two

%{
% Plotting originals
figure, PN_imshowlabel(p,frame1recentered,[],[],[]);
figure, PN_imshowlabel(p,frame2,[],[],[]);
%}

frame1filtered = frame1recentered .* (frame2>0);
%frame2filtered = (frame1recentered>0) .* frame2;

if debugmode
    p.showNr = 1;
    figure, PN_imshowlabel(p,frame1filtered,[],[],[]);
    %figure, PN_imshowlabel(p,frame2filtered,[],[],[]);
end

%% Retrieving the links 

% loop over all the cells, use indices from frame 2 as reference
linklist = [];
for cellidxsin2 = 1:max(frame2(:)) % by construction these are indices 
% Note that "1" and "2" suffixes relate to frame 1 and frame 2 resp.
    % retrieve overlap info
    currentcellidx2 = find(frame2==cellidxsin2);
    %frame1filtered(currentcell)
    %unique(frame1filtered(currentcell))
    
    % Quantify the overlap w. respect to previous frame.
    [pixcounts,possibleParentsList] = hist(frame1filtered(currentcellidx2),unique(frame1filtered(currentcellidx2)));
    
    % Get the most likely one
    [maxpixcount, temp_idx] = max(pixcounts);
    mostlikelyparentidx1 = possibleParentsList(temp_idx);
    
    %disp(['Found link ' num2str(cellidxsin2) ' to ' num2str(mostlikelyparentidx1) '.']);
    linklist = [linklist; mostlikelyparentidx1,cellidxsin2]; % [parent, daughter; ..]
end

if debugmode
    linklist

    figure, PN_imshowlabel(p,frame1recentered,[],[],[]);
    figure, PN_imshowlabel(p,frame2,[],[],[]);
end
    
%% Make an origin picture    
close all;

parentlist=linklist(:,1)
daughterlist=linklist(:,2)
frame2parents = changem(frame2, parentlist, daughterlist );

if debugmode
    figure, PN_imshowlabel(p,frame1recentered,[],[],[],'CustomColors',myColors);
    figure, PN_imshowlabel(p,frame2parents,[],[],[],'CustomColors',myColors);
end


%% Convert to schnitzcells format 
    
linklistschnitz = [linklist(:,1), zeros(size(linklist,1),2), linklist(:,2)];
[parentcount,indices] = hist(linklistschnitz(:,1),unique(linklistschnitz(:,1)));

parentlist = indices(find(parentcount>1))';
cellListFrame1 = linklist(:,1);
dividesInLinkList = find(changem(cellListFrame1,parentcount,indices)>1);

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
        disp('ERROR: parent has more than 2 daughters! Leaving untouched!');
    end
end
    
if debugmode
    dividingrowsafter = linklistschnitz(dividesInLinkList,:)    
end
    
    
%% Perform check(s)

checksPassed = 1;

orphans = find(linklist(:,1)==0);
if ~isempty(orphans)
    disp('ERROR: orphan cells found. At lines:');
    orphans
    checksPassed = 0;
end

barren = find(linklist(:,2)==0);
if ~isempty(barren)
    disp('ERROR: barren cells found. At lines:');
    barren
    checksPassed = 0;
end

% Daughters of wich ancestry is contended (should be impossible)
% This is impossible by construction (hence commented out).
%{
[daughtercount,indices] = hist(linklistschnitz(:,4),unique(linklistschnitz(:,4)));
milkmanKids = indices(find(daughtercount>1));
if ~isempty(milkmanKids)
    disp('Contented ancestry found.')
    milkmanKids
    checksPassed = 0;
end
%}

if checksPassed 
    disp('All checks passed..')
end
    
    
    
    
    
    
    

