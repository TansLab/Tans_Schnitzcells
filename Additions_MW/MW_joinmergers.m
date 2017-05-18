function  affectedcells = MW_joinmergers(p,L,rect,Lp,rectp)
% During segementation correction this function allows you to join two
% cells that are detected as one cell in the next frame.
%
% input 
% p     standard input
% L     segmentation image
% rect  crop recteangle of segm image
% Ln    segmentation image of next frame
% rectn crop rect of next frame
%
% optional input
% p.extrashow
%
% This function is totally not optimized and might do things that are just
% wasting CPU time. (TODO)

%% Load some test images
%{
load('H:\EXPERIMENTAL_DATA_2017\2017-03-22_asc1004_cAMP_pulsing\pos2smallcrop\segmentation\pos2smallcropseg1369.mat');
Lp=Lc; rectp=rect;
load('H:\EXPERIMENTAL_DATA_2017\2017-03-22_asc1004_cAMP_pulsing\pos2smallcrop\segmentation\pos2smallcropseg1370.mat');
L=L; rect=rect;
%}

%% code below is stolen from PN_imshowlabel

% get relative movement using average colony centers in the 2 frames
[rw clm] = find(Lp);
center_old = round([mean(clm) mean(rw)]);
[rw clm] = find(L);
center_new = round([mean(clm) mean(rw)]);
pad_motion = center_new + [rect(2) rect(1)] - [rectp(2) rectp(1)] - center_old; 

% get the centroids from previous image
propLp = regionprops(Lp,'Centroid');
centroidsraw = cat(1, propLp.Centroid);
centroids = round(centroidsraw + ones(size(propLp)) * ([rectp(2)-1 rectp(1)-1] + pad_motion)); % add offset from cropping and path motion
imcentroids = zeros(max(rect(3),rectp(3)),max(rect(4),rectp(4)));
    
try
    linearInd = sub2ind(size(imcentroids), centroids(:,2), centroids(:,1));
catch
    warning('Converting centroids to linear indices failed.');
    linearInd = NaN;
end

imcentroids(linearInd) = 1;

imcentroids = imcentroids(rect(1):rect(3),rect(2):rect(4));
%find domains which have 2 or more former centroids
Ltemp = L;
Ltemp(logical(imcentroids)) = 0;
propLtemp = regionprops(Ltemp,'EulerNumber');
ideul = find([propLtemp.EulerNumber]<0);

%% 
if isfield(p,'extraShow');
    p.showPhaseImage=0;
    figure(1); PN_imshowlabel(p,Lp,[],[],[]);
    hold on; plot(centroidsraw(:,1),centroidsraw(:,2),'ow');
    figure(2); PN_imshowlabel(p,L,[],[],[]);
    hold on; plot(centroidsraw(:,1),centroidsraw(:,2),'ow');
end

%% look for cells that are merging back to one cell

% for reference create table that gives cellno(centroidnr)
cellnoForCentroid=Lp(sub2ind(size(L),round(centroidsraw(:,2))',round(centroidsraw(:,1))')); % in principle this should be 1:N

% look which cells are hit by centroids
hits = L(sub2ind(size(L),round(centroidsraw(:,2))',round(centroidsraw(:,1))'));
% now count how many times cells got hit by those centroids
duplicity = arrayfun(@(x) sum(x==hits), hits);
% alternative way
%[counts, bins]=histcounts(hits,[-.5:max(L(:))+1])
%duplicity = counts(hits+1)
% now find cells that get hit multiple times
joinedCellIdxs = find(duplicity>1 & hits~=0);
%cellnoForCentroid(joinedCellIdxs)
if isfield(p,'extraShow');
    figure(1); hold on; plot(centroidsraw(joinedCellIdxs,1),centroidsraw(joinedCellIdxs,2),'sw','MarkerFaceColor','w');
end

% now make a list of pairs that need to be connected
buddies = arrayfun(@(x) find(x==hits), hits,'UniformOutput',false);
% ignore cells that are barren
for backgroundHits = find(hits==0)
    buddies{backgroundHits} = [];
end

%% now fix the pairs

Lout=Lp;
%Lout_undo=Lout;

affectedcells=[];
for buddyIdx = 1:numel(buddies)
    
    if numel(buddies{buddyIdx})>1                
        
        buddy1=buddies{buddyIdx}(1);
        buddy2=buddies{buddyIdx}(2);
        
        buddy12=sort([buddy1,buddy2]);
        newColor = buddy12(1);
        oldColor = buddy12(2);
        
        affectedcells=[affectedcells, buddy12]; % just to keep track
        
        pt1=[centroidsraw(buddy1,2),centroidsraw(buddy1,1)];
        pt2=[centroidsraw(buddy2,2),centroidsraw(buddy2,1)];
        
        Lout(Lout==oldColor)=newColor;
        Lout=drawline(Lout, pt1, pt2, newColor);         
        
        buddies{buddy1}=[];
        buddies{buddy2}=[];
        
        disp(['Joined cells ' num2str(buddy1) ' and ' num2str(buddy2) ' in frame ' num2str(p.currentFrame-1) '.']);
    end
    
end

%%
% figure(3); PN_imshowlabel(p,Lout,[],[],[]);

%% Make an additional save of the previous frame
Lc = Lout;
fileLocation=[p.segmentationDir p.movieName 'seg' sprintf('%03d', p.currentFrame-1) '.mat'];
save(fileLocation,'Lc','-append');
disp(['Additional file save to ' fileLocation '.']);













