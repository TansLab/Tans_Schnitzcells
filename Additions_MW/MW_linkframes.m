function [linklistschnitz, segFile1Path, segFile2Path] = MW_linkframes(p, frame1Number, frame2Number)
% function [linklistschnitz, segFile1Path, segFile2Path] = MW_linkframes(p, frame1Number, frame2Number)
%
% Performs tracking for frames frame1Number and frame2Number.
%
% if p.debugmode is valid field, then also figures are plotted.
%
% Note term parent and daughter are here also used to link the same
% individual over two frames!
%
% Returns:
% - segFile1Path, segFile2Path:
%       Location of .mat segmentation files that were linked.
% - linklistschnitz
%       Simple list, where each row links segmentation nr of cell in frame
%       n to segmentation nr of cell in frame n+1.
%       EXCEPTIONS, returns:
%       - linklistschnitz = 0 if frame pair was skipped because seg file
%         older than track file.
%       - linklistschnitz = -1 if checks were not passed.


%% Parameters 

DISKSIZE=15;
MARGIN=10;

if ~exist('frame1Number','var')
    frame1Number=151;
    frame2Number=152;
end

if isfield(p,'debugmode')
    debugmode=p.debugmode;
else
    debugmode=0;
end

if ~isfield(p,'overwrite')
    p.overwrite=0;
end
if isfield(p,'override') % backwards compatibility
    p.overwrite=p.override;
end

% Data file names
myFileStringStart = [p.dateDir p.movieName '\segmentation\' p.movieName 'seg'];
% Filename 1
segFile1Path = [myFileStringStart  sprintf('%03d', frame1Number) '.mat'];
% Filename 2
segFile2Path = [myFileStringStart sprintf('%03d', frame2Number) '.mat'];
% Output file
trackOutputFile = [p.tracksDir,p.movieName,'-djk-output-',str3(frame1Number),'-to-',str3(frame2Number),'.txt'];

%% Quit analysis if tracking file is newer than segfile
% (Code stolen from NW_tracker_centroid_vs_area.)
if exist(trackOutputFile,'file')==2 & ~p.overwrite

    % if trackOutputFile is younger than segmentation files, means that segmentation has been updated, and we want to track again
    info_yesterdaySegFile = dir(segFile1Path);
    info_todaySegFile = dir(segFile2Path);
    info_trackOutputFile = dir(trackOutputFile);
    datenumber_yesterdaySegFile = datenum(info_yesterdaySegFile.date);
    datenumber_todaySegFile = datenum(info_todaySegFile.date);
    datenumber_trackOutputFile = datenum(info_trackOutputFile.date);
    if ~((datenumber_yesterdaySegFile>datenumber_trackOutputFile) || (datenumber_todaySegFile>datenumber_trackOutputFile))
        fprintf(1,' -> Skipping, cause seg older than previous tracking (use p.overwrite=1 to redo)\n');
        linklistschnitz = 0;
        return
    end
end

%% Loading files
% Load file 1
data=load(segFile1Path,'Lc');
frame1=data.Lc;
if debugmode
    figure(1)
    PN_imshowlabel(p,frame1,[],[],[]);
end
% Load file 2
data=load(segFile2Path,'Lc');
frame2=data.Lc;
if debugmode
    figure(2), PN_imshowlabel(p,frame2,[],[],[]);
end

myColors = [0,0,0; distinguishable_colors(max([frame1(:); frame2(:)]),[0,0,0])];

%%
% DISKSIZE defined at top

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
if debugmode
    % Find 
    perimImg = bwperim(frame2BWarea);
    %perimImg = imdilate(perimImg,strel('disk',3,4)); % imdilate use TODO can be optimized

    areaImg1 = frame1BWarea+frame1;
    figure(4), PN_imshowlabel(p,areaImg1,[],[],[]);

    areaImg2 = frame2BWarea+frame2;
    figure(5), PN_imshowlabel(p,areaImg2,[],[],[]);    
    
end

%%
STATS1 = regionprops(frame1BWarea, 'Area');
STATS2 = regionprops(frame2BWarea, 'Area');
    
if (numel(STATS1) > 1) || (numel(STATS2) > 1)
    warning('Failed to detect single colony area, attempting to fix by summing areas.');
    sum1 = sum(STATS1.Area);
    sum2 = sum(STATS2.Area);
    clear STATS1; STATS1.Area = sum1;
    clear STATS2; STATS2.Area = sum2;
end

resizeratio = sqrt(STATS2.Area/STATS1.Area);
    
frame1resized = imresize(frame1, resizeratio,'nearest');

if debugmode
    resizeratio
    figure(6), PN_imshowlabel(p,frame1resized,[],[],[]);
end

 
%% create aligned frame #1
% MARGIN defined at top

% Create dilate/erode frame1resized
frame1resizedBW = (frame1resized>0);
frame1resizedBWarea=imdilate(frame1resizedBW,se);
frame1resizedBWarea=imerode(frame1resizedBWarea,se);
% Find centroids
STATS1 = regionprops(frame1resizedBWarea, 'Centroid');
STATS2 = regionprops(frame2BWarea, 'Centroid');
if (numel(STATS1) > 1) || (numel(STATS2) > 1)
    warning('Failed to detect single colony area. Attempting to fix by averaging centroids.');
    meanx1 = mean(arrayfun(@(k) STATS1(k).Centroid(1), 1:numel(STATS1)));
    meany1 = mean(arrayfun(@(k) STATS1(k).Centroid(2), 1:numel(STATS1)));
    meanx2 = mean(arrayfun(@(k) STATS2(k).Centroid(1), 1:numel(STATS2)));
    meany2 = mean(arrayfun(@(k) STATS2(k).Centroid(2), 1:numel(STATS2)));   
    clear STATS1; STATS1.Centroid = [meanx1, meany1];
    clear STATS2; STATS2.Centroid = [meanx2, meany2];
end

% Create frame that has equal size to frame 2, but contains info from
% frame1, with centroid at same position as frame 2.

sizeframe1 = size(frame1resized);
sizeframe2 = size(frame2);

deltaSize = sizeframe2-sizeframe1;
%deltaSize(find(deltaSize<0)) = 0; % if shrinks, shouldn't happen

frame1recentered = zeros(sizeframe2);

% determine translation
deltacentroids = round(STATS2.Centroid-STATS1.Centroid);
%deltacentroids = [0,0]; % MW DEBUG REMOVE
deltacentroidi = deltacentroids(2); deltacentroidj = deltacentroids(1);

% Get the coordinates for the resized frame 1
iInEnlarged = [1,sizeframe1(1)]+deltacentroidi; % MW grr
jInEnlarged = [1,sizeframe1(2)]+deltacentroidj; % MW grr
iInOriginal = [1,sizeframe1(1)]; 
jInOriginal = [1,sizeframe1(2)];

% Some horrible administration for when translated falls outside window
if iInEnlarged(1)<1, iInEnlarged(1)= 1; iInOriginal(1)= 1-deltacentroidi;  end
if iInEnlarged(2)>sizeframe2(1), iInEnlarged(2)= sizeframe2(1); iInOriginal(2)= sizeframe2(1)-deltacentroidi;  end
if jInEnlarged(1)<1, jInEnlarged(1)= 1; jInOriginal(1)= 1-deltacentroidj;  end
if jInEnlarged(2)>sizeframe2(2), jInEnlarged(2)= sizeframe2(2); jInOriginal(2)= sizeframe2(2)-deltacentroidj;  end

if debugmode
    sizeframe1
    sizeframe2
    iInEnlarged, iInOriginal
    jInEnlarged, jInOriginal
    deltacentroids
end    

frame1recentered(iInEnlarged(1):iInEnlarged(2), jInEnlarged(1):jInEnlarged(2)) = ...
    frame1resized(iInOriginal(1):iInOriginal(2), jInOriginal(1):jInOriginal(2));

% Recenter image
%{
frame1recentered = padarray(frame1,deltaSize,'post'); % resize to have equal size
frame1recentered = imtranslate(frame1recentered,deltacentroids); % align
%}
             
             
%% plot overlay of aligned images
             
if debugmode
    frameBW = (frame1recentered>0);
    frameBWarea=imdilate(frameBW,se);
    frameBWarea=imerode(frameBWarea,se);
    
    overlayareas = frameBWarea+frame2BWarea;
    figure, PN_imshowlabel(p,overlayareas,[],[],[]);
    
    overlaycentered = frame1recentered+frame2;
    figure, PN_imshowlabel(p,overlaycentered,[],[],[]);
    
    
    % plotting areas
    m1=(frameBWarea>0)*1;
    m2=(frame2BWarea>0)*2;
    m=m1+m2;
    figure, PN_imshowlabel(p,m,[],[],[]);
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
    uniqueExclZeros = unique(frame1filtered(currentcellidx2));
    uniqueExclZeros = uniqueExclZeros(find(uniqueExclZeros>0));   
    if numel(uniqueExclZeros)>1
        [pixcounts,possibleParentsList] = hist(frame1filtered(currentcellidx2),uniqueExclZeros);
    
        % Get the most likely one
        [maxpixcount, temp_idx] = max(pixcounts);
        mostlikelyparentidx1 = possibleParentsList(temp_idx);
    else
        mostlikelyparentidx1 = uniqueExclZeros;
        if isempty(mostlikelyparentidx1)
            mostlikelyparentidx1 = 0;
        end
    end
    
    %disp(['Found link ' num2str(cellidxsin2) ' to ' num2str(mostlikelyparentidx1) '.']);
    linklist = [linklist; mostlikelyparentidx1,cellidxsin2]; % [parent, daughter; ..]
end

if debugmode
    linklist

    %figure, PN_imshowlabel(p,frame1recentered,[],[],[]);
    %figure, PN_imshowlabel(p,frame2,[],[],[]);
    
    disp('Section done');
end
    
%% Make an origin picture    
parentlist=linklist(:,1);
daughterlist=linklist(:,2);
frame2parents = changem(frame2, parentlist, daughterlist );

if debugmode
    figure, PN_imshowlabel(p,frame1recentered,[],[],[],'CustomColors',myColors);
    figure, PN_imshowlabel(p,frame2parents,[],[],[],'CustomColors',myColors);
end

%% Perform check(s)

checksPassed = 1;

orphans = find(linklist(:,1)==0);
if ~isempty(orphans)
    disp('ERROR: orphan cells found. At lines:');
    orphansinlinklistindex=orphans
    checksPassed = 0;
end

%barren = find(linklist(:,2)==0);
uniqueExclZerosFr1 = unique(frame1); % get list w. cellno's in this frame
uniqueExclZerosFr1 = uniqueExclZerosFr1(find(uniqueExclZerosFr1>0)); % but idx=0 isn't a cell, filter out 
Frame1LinkedOnes = ismember(uniqueExclZerosFr1, linklist(:,1)); % check whether the cellno's from 1 are all accounted for in list
barren = find(Frame1LinkedOnes==0); % if not, they're barren
if ~isempty(barren)
    disp('WARNING: barren cells found. cellno''s:');
    barren
    checksPassed = 0;
    disp('This leads to serious issue, since leads to discrepancy between numel(schnitzcells(i).frame_nrs) and numel(schnitzcells(i).cellno)!');
    %warning('Adding them as connected to schnitz #1..');
    %linklist = [linklist; padarray(barren, [0,1],1,'post')];
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
else
    disp(['WARNING: Checks not passed. Skipping frames ' num2str(frame1Number) '-' num2str(frame2Number) '. Re-track with other tracker.']);
        % other trackers: DJK_tracker_djk or NW_tracker_centroid_vs_area
    linklistschnitz = -1;
    return
end

%% Convert to schnitzcells format 
   
linklistschnitz=MW_convertsimplelinklisttoschnitzlist(linklist,debugmode); 

    
    
%% Writing output

% Code re. trackOutputFile stolen from: NW_tracker_centroid_vs_area
% ===
%{
if (exist(trackOutputFile) == 2) && (p.overwrite == 0)
    disp('ERROR: didn''t write output since file already exist (and p.overwrite == 0)');
else
%}

% clean up any existing output
if exist(trackOutputFile) == 2
    disp('WARNING: deleted old trackingfile.');
    delete(trackOutputFile)
end

% Open trackOutputFile
fid = fopen(trackOutputFile,'wt');

% loop over results and print to file
for i = 1:length(linklistschnitz(:,1))
    fprintf(fid, '%u %u %u %u\n', linklistschnitz(i,1), linklistschnitz(i,2), linklistschnitz(i,3), linklistschnitz(i,4)); 
end

% Close trackOutputFile
fclose(fid);
    
%{    
end
%}
    
    
disp('Finished');    
    











