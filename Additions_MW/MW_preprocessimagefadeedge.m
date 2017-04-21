
function theImage = MW_preprocessimagefadeedge(theImage, bottomOrTop)
% This function adds a band of maximum intensity at the edge of an image,
% this band has a certain thickness, over which it fades away
%
% This is an addition that could be considered somewhat of a
% "hack". It pre-processes images by adding a gradient at one edge
% to make sure cells that are at the edge of the image are
% processed conveniently. 
% This leads to cells that are partially detected at the edges,
% which are removed later by MW_deletecellsattheedge.
%
% TODO: Note that currently the band is added at the bottom of the image,
% this should be changed such that it can also be added to the top with a
% simple option..

%%
%testImg2=imread('G:\EXPERIMENTAL_DATA_2017\2017-02-10_OAA_pulsing_asc1004\pos2crop\images\pos2crop-p-2-2678.tif');

maxValue = max(theImage(:));

% fade the bottom
if bottomOrTop==1
    theEnd = size(theImage,1);


    THICKNESS = 20;
    for i=0:THICKNESS

        theImage(theEnd-i, :) = (1-i/THICKNESS) .* double(maxValue) .* ones(1,size(theImage,2)) +....
                                (  i/THICKNESS) .* double(theImage(theEnd-i, :));
    end
% or the top    
elseif bottomOrTop==2
    THICKNESS = 20;
    for i=0:THICKNESS

        theImage(i+1, :) = (1-i/THICKNESS) .* double(maxValue) .* ones(1,size(theImage,2)) +....
                           (  i/THICKNESS) .* double(theImage(i+1, :));
    end
    
else
    error('Unrecognized case of bottomOrTop.');
end



%{
figure; imshow(theImage,[]);
imageToSegment = theImage;
%}
