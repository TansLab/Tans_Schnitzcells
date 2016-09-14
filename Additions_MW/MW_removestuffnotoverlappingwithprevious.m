
function [Lout, Lout_undo] = MW_removestuffnotoverlappingwithprevious(p,Lout,L_prec,rect,rect_prec)
%
% This function removes stuff that does not overlap with cells as
% from previously corrected image.
%
% Arguments:
% Lout          input image
% L_prec        image from preceeding frame
% rect          ROI recteangle
% rect_prec     rect from preceeding img
%
% Original code written by Noreen Walker. Edited by Martijn Wehrens.

DILATIONPERCENTAGE = .01;
DEBUGMODE = 0;

%% For debug purposes
if DEBUGMODE
    MYFRAMENR1=15;  % "previous" frame
    MYFRAMENR2=16;  % "current" frame
    % load params previous frame
    filename = [p.segmentationDir,p.movieName,'seg',str3(MYFRAMENR1),'.mat']; 
    load(filename);
    % rename for testing
    L_prec = Lc;
    rect_prec = rect;
    % load params
    filename = [p.segmentationDir,p.movieName,'seg',str3(MYFRAMENR2),'.mat']; 
    load(filename);
    % rename for testing
    Lout = LNsub;
    % other params needed
    p.fullsize = [1024 , 1344];
    warning('Test mode engaged.');
end

%% for undo step NW 
Lout_undo=Lout;   

%% old and new seg need to be located properly:
Lout_full=zeros(p.fullsize);
Lout_full(rect(1):rect(3),rect(2):rect(4))=Lout;
L_prec_full=zeros(p.fullsize);
L_prec_full(rect_prec(1):rect_prec(3),rect_prec(2):rect_prec(4))=L_prec;

%% create binary mask
L_prec_mask=L_prec_full;             
%L_prec_mask=imdilate(L_prec_full,strel('disk',3)); % enlargen prev. image overlap mask if needed
L_prec_mask=L_prec_mask>0;  % binary mask from previous image

%% Determine cells to delete

% dilate colony by +/- 5%
% ===
% determine parameters
colonyWidth = sqrt(sum(L_prec_mask(:)));
dilateMaskSize = round(colonyWidth*DILATIONPERCENTAGE);
SE = strel('disk',dilateMaskSize);
% dilate
L_prec_mask_dilated = imdilate(L_prec_mask,SE);               

% multiply previous w. current, such that
% overlapping stuff remains
OverlappingWithPrevious = Lout_full.*L_prec_mask_dilated;

indicesAllCells=unique(Lout)';
indicesOverlappingWithPrevious = unique(OverlappingWithPrevious)';
cellstodelete=setdiff(indicesAllCells,indicesOverlappingWithPrevious);

if DEBUGMODE
    
    % colored version    
    frame1Pic = (1-(Lout_full>0)*.5);
    frame2Pic = (1-(L_prec_mask_dilated>0)*.5);
    frame1PicColored = ones([size(frame1Pic),3]);
    frame1PicColored(:,:,1) = frame1Pic; % color blue channel
    frame1PicColored(:,:,2) = frame1Pic; % color blue channel
    frame2PicColored = ones([size(frame2Pic),3]);
    frame2PicColored(:,:,2) = frame2Pic; % color red channel
    frame2PicColored(:,:,3) = frame2Pic; % color red channel
    % show them
    h1=figure(); imshow(frame1PicColored.*frame2PicColored)
    
    %h2=figure(); imshow(OverlappingWithPrevious,[]);
    %h3=figure(); imshow(Lout,[]);
    cellstodelete
end

%% Actually remove all "cells" that are not (completely) overlapping
Lout(ismember(Lout,cellstodelete)) = 0;

if DEBUGMODE
    h4=figure(); imshow(Lout,[]);
end

