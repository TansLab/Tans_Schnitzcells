
% ===
% THIS FUNCTION IS DEPRICATED AND SHOULD NOT BE USED ANY MORE
% ===
%
% Use PN_imshowlabel instead. That fn has been modified such that it should
% be equivalent to DJK_imshowlabel when its vars rect, Lp and rectp are set 
% to 0.
%
% I.e. use 
% PN_imshowlabel(p,L,[],[],[])
% instead of
% DJK_imshowlabel(L) or DJK_imshowlabel(p,L)
%
% MW, 2015/01



































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
%

function outim = DJK_imshowlabel(p,L,varargin);

disp('THIS FUNCTION (DJK_imshowlabel) IS DEPRICATED AND SHOULD NOT BE USED ANY MORE, USE PN_imshowlabel INSTEAD.')

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'DJK_imshowlabel'; p_internal = struct;

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)))
  errorMessage = sprintf('%s\n%s',['Error with input arguments of ' functionName],['Try "help ' functionName '".']);
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
% Overwrite any schnitzcells parameters/defaults given optional fields/values
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

if addPhaseImage

    rgb = 0.5 * ind2rgb(L2,mymap);
    bwscreen = double(p_internal.phaseImage); % bwscreen = 0.5 * bwscreen / max(max(bwscreen));
    bwscreen = DJK_scaleRange(bwscreen, [max(max(bwscreen)) min(min(bwscreen))], [0 1]);
    bwscreen = DJK_scaleRange(bwscreen, [0.25 1], [0 0.5]);

    rgb(:,:,1) = rgb(:,:,1) + bwscreen;
    rgb(:,:,2) = rgb(:,:,2) + bwscreen;
    rgb(:,:,3) = rgb(:,:,3) + bwscreen;

    outim = rgb;  
    outim = MW_stampit(outim,p); % Edit MW - adds green marker if frame has been approved   
    
    if nargout == 0
        imshow(outim);
    end

else
    
    outim = L2;
    outim = MW_stampit(outim,p); % Edit MW - adds green marker if frame has been approved 

    if nargout == 0
        imshow(outim, mymap);
    end  

end


     
