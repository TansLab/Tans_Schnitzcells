




% Experimental function! 
% Do not use!

% MW 2015/07















FRAMENUMBER = 190;

%function connections = MW_tracker(frame1, frame2)

figure(1)
data=load(['F:\A_Tans1_step1_incoming_not_backed_up\2015-06-02\pos7crop\segmentation\pos7cropseg'  sprintf('%03d', FRAMENUMBER) '.mat'],'LNsub');
frame1=data.LNsub;
PN_imshowlabel(p,frame1,[],[],[]);

figure(2)
data=load(['F:\A_Tans1_step1_incoming_not_backed_up\2015-06-02\pos7crop\segmentation\pos7cropseg' sprintf('%03d', FRAMENUMBER+1) '.mat'],'LNsub');
frame2=data.LNsub;
PN_imshowlabel(p,frame2,[],[],[]);

%%
DISKSIZE=15;

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

figure(4)
PN_imshowlabel(p,areaImg,[],[],[]);

%%
STATS1 = regionprops(frame1BWarea, 'Area', 'Centroid', 'EulerNumber');
STATS2 = regionprops(frame2BWarea, 'Area', 'BoundingBox', 'Centroid', 'EulerNumber');
    
if (numel(STATS1) > 1) || (numel(STATS2) > 1)
    error('Failed to detect single colony area.');
end

ratio = STATS2.Area/STATS1.Area
    
frame1resized = imresize(frame1, ratio);
figure(6), imshow(frame1resized)
   
%% crop
MARGIN=10;

% Create frame that has equal size to frame 2, but contains info from
% frame1, with centroid at same position as frame 2.

frame1recentered = zeros(size(frame2));
deltacentroids = round(STATS2.Centroid-STATS1.Centroid);
sizeframe1 = size(frame1);
frame1recentered([1:sizeframe1(1)]+deltacentroids(1), ...
                 [1:sizeframe1(2)]+deltacentroids(2)) = frame1;

%%
             
overlaycentered = frame1recentered+frame2;
figure, PN_imshowlabel(p,overlaycentered,[],[],[]);

%%
    
    

