
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
%theImage=imread('H:\EXPERIMENTAL_DATA_2017\2017-03-22_asc1004_cAMP_pulsing\pos2smallcrop\images\pos2smallcrop-p-2-1344.tif'); theImageBefore=theImage;
%theImage=imread('H:\EXPERIMENTAL_DATA_2017\2017-03-22_asc1004_cAMP_pulsing\pos2smallcrop\images\pos2smallcrop-p-2-1540.tif'); theImageBefore=theImage;
%MYFILTER = [.0000 0.000 0.000 0.0025 0.0067 0.0180 0.0474 0.1192 0.2689 0.5000 0.7311 0.8808  0.9526  0.9820 0.9933  1  1  1  1  1];    
    % above is modified version of: 1./(1+(exp(-[-10:10]))); %plot(1./(1+(exp(-[-10:10]))),'o')
    % previously: [0:20]/THICKNESS; %THICKNESS = 20;
MYFILTER = [0 0 0 0 0 linspace(0,1,10)  1  1  1  1  1]; 
    % Note that a few completely white pixels (=0 in array above) are
    % convenient to make sure there is really no touching to the edge.
    % Also some fully transparent pixels are convenient to make sure the
    % cells will eventually touch the "removal box" in 
    % MW_deletecellsattheedge().       

maxValue = max(theImage(:));

% fade the bottom
if bottomOrTop==1
  theEnd = size(theImage,1);


  
  for i=0:numel(MYFILTER)-1 %0:THICKNESS

    a = MYFILTER(i+1); % alpha value
    theImage(theEnd-i, :) = (1-a) .* double(maxValue) .* ones(1,size(theImage,2)) +....
                ( a) .* double(theImage(theEnd-i, :));
  end
% or the top  
elseif bottomOrTop==2
  
  for i=0:numel(MYFILTER)-1 %THICKNESS

    a = MYFILTER(i+1); % alpha value
    theImage(i+1, :) = (1-a) .* double(maxValue) .* ones(1,size(theImage,2)) +....
                        ( a) .* double(theImage(i+1, :));
  end
  
else
  error('Unrecognized case of bottomOrTop.');
end



%{
figure; imshow(theImageBefore,[]);
figure; imshow(theImage,[]);
%}
