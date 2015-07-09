% PN_imshowlabel is a modified version of DJK_imshowlabel which provides
% visual aid for segmentation correction
% - too small cells are shown in white
% - cells with potential overlap with 2 cells in the former image are
% shown in white
%
% DJK_imshowlabel is used to display an integer image
% Uses different kind of color map as before
%
% OUTPUT
% 'outim'             color image
%
% REQUIRED ARGUMENTS:
% 'L'                 seg image
%
% OPTIONAL ARGUMENTS:
% 'phaseImage'        phase image of same size as seg, will be shown as
%                     background
% 'randomize' = 0     no randomizing of colormap (default:1)
% TODO: 'recalcCellNumbers' array for which region properties have to be
%                     recalculated because cell was updated. default: all
%                     cells. taking fewer cells will speed up process


% elapsed time for frame 444 in 2012-05-08. 390 cells


function outim = PN_imshowlabel(L,rect,Lp,rectp,varargin)

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 4; functionName = 'PN_imshowlabel'; p_internal = struct;

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)))
    errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
    error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
            error(errorMessage);
        end
        fieldName = DJK_schnitzfield(varargin{i});
        p_internal.(fieldName) = varargin{i+1};
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
addPhaseImage = false;
if existfield(p_internal,'phaseImage') & length(p_internal.phaseImage)>0
    addPhaseImage = true;
end
if ~existfield(p_internal,'randomize')
    p_internal.randomize = 1;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

if islogical(L),
    L = double(L);
end;

%start time (calculate cost)
%tic

%Suspicious cells detection.
%Detect cells that have certainly been divided in the preceeding image
[rw clm] = find(Lp);
center_old = round([mean(clm) mean(rw)]);
[rw clm] = find(L);
center_new = round([mean(clm) mean(rw)]);
pad_motion = center_new + [rect(2) rect(1)] - [rectp(2) rectp(1)] - center_old; 
%Positions of the centroids of the former image in the new image
propLp = regionprops(Lp,'Centroid');
centroids = cat(1, propLp.Centroid);
centroids = round(centroids + ones(size(propLp))*([rectp(2)-1 rectp(1)-1] + pad_motion));
imcentroids = zeros(max(rect(3),rectp(3)),max(rect(4),rectp(4)));
linearInd = sub2ind(size(imcentroids), centroids(:,2), centroids(:,1));
imcentroids(linearInd) = 1;
imcentroids = imcentroids(rect(1):rect(3),rect(2):rect(4));
%find domains which have 2 or more former centroids
Ltemp = L;
Ltemp(logical(imcentroids)) = 0;
propLtemp = regionprops(Ltemp,'EulerNumber');
ideul = find([propLtemp.EulerNumber]<0); %label of cells which may have been re-merged
%detect too small cells
propL = regionprops(L,'Area');
characSize = median([propL.Area]);
idsmall = find([propL.Area] < characSize/2); %label of cells which may be too small
%create a logical of suspicious cells
Lsuspicious = zeros(size(L));
combinedcellidx = union(ideul,idsmall); % do below loop only when non-empty (NW2015-05. bug appeared with new matlab version?)
if ~isempty(combinedcellidx)
    for ii=combinedcellidx
        Lsuspicious(L==ii)=1;
    end
end

% MW: make them stand out some more using checkerboard
sqsize = 4; %square size
mycheckerboard = (checkerboard(sqsize,ceil(size(Lsuspicious,1)/sqsize),ceil(size(Lsuspicious,2)/sqsize)) > 0.5);
for i = 1:size(Lsuspicious,1)
    for j = 1:size(Lsuspicious,2)
        if mycheckerboard(i,j)==0
            Lsuspicious(i,j) = 0;
        end
    end
end

% elapsed time: 0.33
%stop1=toc

%creation of a rgb colormap
% L2 has every non-background blob in range [2,256] and 
% sets background to one, corresp. to first entry in mymap
L2 = mod(L,255)+2;
L2(L==0) = 1;
% M is the maximum color table entry, at most 256 colors
M = min(max2(L)+2,256);
% create a color map
mymap = DJK_hsv(M); % DJK 071207
% explicitly set the colormap's first entry to black for background
mymap(1,:)=[0 0 0];
if p_internal.randomize
  % get sequence of random integers in range [1,maxcolors-1]
  [s,I] = sort(rand(M-1,1));  
  % randomly reorder mymap color entries [2,maxcolors]
  mymap(2:end,:) = mymap(I+1,:);
end
mymap = [mymap ; 1 1 1]; %add white

% elapsed time: 0.34
%stop2=toc
L2(logical(Lsuspicious)) = M+1;
L2(logical(imdilate(imcentroids,strel('square',2)))) = 0; % costs 0.017 sec
Lrgb = ind2rgb(L2,mymap); % costs 0.04 sec

% elapsed time: 0.40
%stop3=toc
if addPhaseImage % costs 0.045 sec
    rgb = 0.5 * Lrgb;
    bwscreen = double(p_internal.phaseImage); % bwscreen = 0.5 * bwscreen / max(max(bwscreen));
    bwscreen = DJK_scaleRange(bwscreen, [max(max(bwscreen)) min(min(bwscreen))], [0 1]);
    bwscreen = DJK_scaleRange(bwscreen, [0.25 1], [0 0.5]);
    rgb(:,:,1) = rgb(:,:,1) + bwscreen;
    rgb(:,:,2) = rgb(:,:,2) + bwscreen;
    rgb(:,:,3) = rgb(:,:,3) + bwscreen;
    outim = rgb;
    
    % elapsed time: 0.45
    %stop4a=toc
else
    outim = Lrgb;
    
    %stop4b=toc
end

%elapsed time: 0.45
%stop5=toc
