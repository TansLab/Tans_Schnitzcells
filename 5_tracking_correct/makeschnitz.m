function p = makeschnitz(p, varargin)
% MAKESCHNITZ     Make schnitzcells lineage structure from original tracking 
% 
%   MAKESCHNITZ converts an 'original' tracking program output to the 
%   schnitzcells cell lineage representation.  This function is given a 
%   schnitzcells parameter structure p that provides the moveie name, track 
%   directory, and segmentation directory.  The program  will look for an 
%   original tracking output file named [p.trackDir p.movieName '_Tdata.mat'].
%   
%   MAKESCHNITZ(P,'Field1',Value1,'Field2',Value2,...) also converts an 
%   original tracking output to schnitzcells, but the extra arguments permit 
%   users to adjust any parameters describing the movie or parameters 
%   controlling the cell tracking process.  The extra arguments can override 
%   any specific parameters provided in P by setting P.Field1 = Value1, 
%   P.Field2 = Value2, etc.  Thus any/all schnitzcells parameter values can 
%   be defined in the function call via these optional field/value pairs.  
%   (This is in the style of setting MATLAB properties using optional 
%   property/value pairs.)  
%   
%   MAKESCHNITZ produces an output file (a MATLAB binary file) that 
%   contains the 'schnitzcells' variable describing each cell's lineage. The 
%   path name of this file defaults to [p.tracksDir,p.movieName,'_lin.mat'], 
%   but can be specified by providing an optional 'lineageName' field in the 
%   schnitcells parameter structure, or by providing an extra 
%   'lineageName',value pair of arguments to this function. 
%   
%   MAKESCHNITZ returns a struct (1x1 struct array) referred to as the 
%   schnitzcells parameter structure that contains fields and values 
%   contained in p, including unchanged/original parameters plus any of those 
%   added or overridden in the list of properties & values provided in the 
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control makeshnitz
% 
%   trackName     file name of track file, defaults to 
%                 [p.trackDir p.movieName '_Tdata.mat'].
%   
%   trackRange    range of frame numbers to extract from the original tracking 
%                 results and convert to schntizcells; by default, all frames 
%                 contained in the original tracking results will be converted; 
%   
%   lineageName   the schnitzcells cell lineage structure is written to this 
%                 file, by default named [p.trackDir p.moveiName '_lin.mat'].
%   
%
%-------------------------------------------------------------------------------
% Varibles contained in the output file:
%
%  schnitzcells   schnitzcells is a 1xNS struct array (NS = number of "schnitz" 
%                 tracks).  This struct array describes the lineage of all 
%                 cells in a movie.  A "schnitz" is a track that begins at 
%                 the first frame a distint cell appears in, and ends at the 
%                 last frame that distinct cell appears in.  When a cell 
%                 divides it's track ends, and two new tracks being at the 
%                 first frame it's two children appear in.  A schnitz (track) 
%                 is numbered by its position (index) within the struct array. 
%                 A minimal schnitz structure is described by P (the schnitz 
%                 number of it's parent), children schnitz numbers D and E,
%                 sister schnitz number S, frames (an array) listing each frame 
%                 this schnitz occurs in, cellno (an array) listing each 
%                 cell number of the tracked cell within each of the frames, 
%                 and N, the number of frames that distinct cell appears.
%                 
%-------------------------------------------------------------------------------% Example inputs:
%   
%   p.movieName = 'x1TimeLapse-11';
%   p.trackRange = [0:115];
%   p.segmentationDir = 'G:\movie_analysis\2005-03-03\x1TimeLapse-11\segmentation\';
%   p.trackName='x1TimeLapse-11_Tdata.mat';
%   p.lineageName = 'G:\movie_analysis\2005-03-03\x1TimeLapse-11\data\x1TimeLapse-11_lin.mat';


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
   (mod(nargin,2) == 0) | ...
   (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
      'Error using ==> makeschnitz:',...
      '    Invalid input arguments.',...
      '    Try "help makeschnitz".');
    error(errorMessage);
end

%-------------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values 
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
              'Error using ==> makeschnitz:',...
              ' Invalid property ', num2str(varargin{i}), ...
              ' is not (needs to be) a string.',...
              '    Try "help makeschnitz".');
            error(errorMessage);
        end
        fieldName = schnitzfield(varargin{i});
        p.(fieldName) = varargin{i+1};
    end
end


if ~existfield(p,'trackName')
    p.trackName = [p.tracksDir,p.movieName,'_Tdata.mat'];
end

if ~existfield(p,'lineageName')
    p.lineageName = [p.tracksDir,p.movieName,'_lin.mat'];
end

%-------------------------------------------------------------------------------
% create schnitzcells lineage structure given tracking results
%-------------------------------------------------------------------------------

% Get track matrix CellMat
trackdata = load(p.trackName);
CellMat = trackdata.cellmat;

numFrames = size(CellMat,1);
currS = 1;

if existfield(p,'trackRange')
    % let user-specified track range select subset
    trackRange = p.trackRange;  
    cellmatrows = ismember(trackdata.trackRange,p.trackRange);
    CellMat = CellMat(cellmatrows,:);
    numFrames = size(CellMat,1);
else
    trackRange = trackdata.trackRange;
end


% Get coordinates of cells in first image
segdata = load([p.segmentationDir,p.movieName,'seg',str3(trackRange(1))]);
pts = regionprops(segdata.Lc,'Centroid');

% Initialize schnitcells from first row
previous = CellMat(1,:);
[temp] = unique(CellMat(1,:));
[prIds] = temp(find(temp ~=0));
SchnitzId = [];
for i = 1:size(prIds,2)
    % Make new Schnitz
    sList(currS).P = -1;
    sList(currS).frames(1) = trackRange(1) + 1; % JCR: add one (for now) to match schnitzedit
    sList(currS).cellno(1) = prIds(i);
    SchnitzId(prIds(i)) = currS;
    % Read in next two from *.mat files if needed
    sList(currS).cenx(1) = pts(prIds(i)).Centroid(1);
    sList(currS).ceny(1) = pts(prIds(i)).Centroid(2);
    sList(currS).D = -1;
    sList(currS).E = -1;
    sList(currS).N = 1;
    currS = currS+1;
end;

for i = 2:numFrames

    % Get coordinates of cells in i imag
    segdata = load([p.segmentationDir,p.movieName,'seg',str3(trackRange(i))]);
    clear pts;
    pts = regionprops(segdata.Lc,'Centroid');

    current = CellMat(i,:);
    clear prevS;
    prevS = SchnitzId(prIds);  %!!! What if it does not exist. Make new schnitz!
    SchnitzId = [];
    for j = 1:size(prIds,2)
        clear prevIn;
        [prevIn] = find(previous == prIds(j));
        clear curIds; 
        curIds = unique(current(prevIn));

        if (size(curIds,2)==1)
            if (curIds(1)~=0)
                % single daughter found
                % Continue j's Schnitz
                sList(prevS(j)).N = sList(prevS(j)).N+1;
                np = sList(prevS(j)).N;
                sList(prevS(j)).frames(np) = trackRange(i) + 1;  % JCR: hack
                sList(prevS(j)).cellno(np) = curIds(1);
                SchnitzId(curIds(1)) = prevS(j);
                % Read in next two from *.mat files if needed
                sList(prevS(j)).cenx(np) = pts(curIds(1)).Centroid(1);
                sList(prevS(j)).ceny(np) = pts(curIds(1)).Centroid(2);
            end;
        else 
            if (size(curIds,2)>=2)
                % (at least) two daughters found
                % set daughters of parent
                sList(prevS(j)).D = currS;
                sList(prevS(j)).E = currS+1;

                % Make two new Schnitz for daughters
                for k = 1:2
                    if (curIds(k)~=0)
                        sList(currS).P = prevS(j);
                        sList(currS).frames(1) = trackRange(i) + 1; % JCR hack
                        sList(currS).cellno(1) = curIds(k);
                        SchnitzId(curIds(k)) = currS;
                        % Read in next two from *.mat files if needed
                        sList(currS).cenx(1) = pts(curIds(k)).Centroid(1);
                        sList(currS).ceny(1) = pts(curIds(k)).Centroid(2);
                        sList(currS).D = -1;
                        sList(currS).E = -1;
                        sList(currS).N = 1;
                        currS = currS+1;
                    end;
                end;

                if (size(curIds,2)>2)
                    % many daughters found
                    sprintf(['Error in frame %d, for previous id %d, '...
                              'illegal number of siblings: %d\n',...
                              i, prIds(j), size(curIds,2)]);
                    % Make new Schnitz for illegal daughters
                    for k = 3:(size(curIds,2))
                        sList(currS).P = -1;
                        sList(currS).frames(1) = i;
                        sList(currS).cellno(1) = curIds(k);
                        SchnitzId(curIds(k)) = currS;
                        sList(currS).D = -1;
                        sList(currS).E = -1;
                        sList(currS).N = 1;
                        currS = currS+1;             
                    end;
                end;
            end;
        end;
    end;

    previous = CellMat(i,:);
    clear prIds;
    [temp] = unique(CellMat(i,:));
    [prIds] = temp(find(temp ~=0)); 
end;

schnitzcells = sList;

disp(['saving schnitzcells lineage structure to ' p.lineageName]);
save(p.lineageName,'schnitzcells');

return;
