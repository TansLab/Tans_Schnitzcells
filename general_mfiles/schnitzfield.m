function standardFieldString = schnitzfield(fieldNameString)
%SCHNITZFIELD Converts Schnitzcells case-insensitive fieldnameString to standard
%
%   SCHNITZFIELD(fieldnameString) returns the schnitzells-standard 
%   "camel-case" schnitzcells structure field name string given the 
%   user-provided fieldNameString that may not have been properly 
%   capitalized.  This utility can be used by any of the schnitzcells 
%   programs to convert 'tainted' user inputs to the standard form before 
%   inserting or updating fields in the schnitzcells parameter structure.
%
%   e.g. schntizfield('rootdir') returns 'rootDir'
%   
%   If the given fieldNameString is not known by this routine, a warning 
%   will be printed, but and the given fieldNameString will be returned. 
%   
%   If a new schntizcells program needs a new field in the schnitzcells 
%   parameter structure, it is recommended practice to add that field 
%   to this function so that a standard capitalization can be recognized.

% Note: it would be nice to one day make this a hash table so that we 
% don't need to do linear search of all fields when matching the given 
% fieldNameString. Hey, for that matter, it would be nice to make it a 
% static hash table (loaded from binary? binary updated if out of date?) 
% rather than reconstructing the fields & upfields info every time it's 
% called. Hmmm... Python...

fields = {...
  'Dskip' ...                       % old manualcheckseg param
  'angThresh' ...
  'autoCFL' ...                     % compileschnitz param
  'autoYFL' ...                     % compileschnitz param
  'autoGFL' ...                     % compileschnitz param
  'autoRFL' ...                     % compileschnitz param
  'beadCFP' ...                     % compileschnitz param
  'beadYFP' ...                     % compileschnitz param
  'beadGFP' ...                     % compileschnitz param
  'beadRFP' ...                     % compileschnitz param
  'cback0' ...                      % compileschnitz param
  'cell' ...                        % fillinator,fill2seg param
  'dateDir' ...
  'edge_lapofgauss_sigma' ...
  'expandvalue' ...                 % old manualcheckseg param
  'fillinatorFile' ...              % fillinator and fill2seg param
  'fillinatorRange' ...             % fillinator param
  'finetuneimage' ...               % old manualcheckseg param
  'flatFieldFile' ...               % compileschnitz param
  'frnum' ...                       % old manualcheckseg param
  'gback' ...                       % compileschnitz param
  'imNumber1' ...
  'imNumber2' ...
  'imageDir' ...
  'leftend' ...                     % old manualcheckseg param
  'len' ...                         % trackcomplete param
  'lineageName' ...                 % schnitzedit param
  'magic' ...                       % compileschnitz param
  'manualRange' ...                 % manualcheckseg param
  'maxCellWidth' ...
  'maxThresh' ...
  'maxThreshCut' ...
  'maxThreshCut2' ...
  'crosstalkRatioC2Y', ...          % compileschnitz param
  'micronsPerPixel', ...            % compileschnitz param
  'minCellArea' ...
  'minCellLength' ...
  'minCellLengthConservative' ...
  'minNumEdgePixels' ...            % segmoviephase param
  'miniMovieName' ...               % makeminimovie param
  'miniRange' ...                   % makeminimovie param
  'minThresh' ...
  'min_size' ...                    % old manualcheckseg param
  'movieDate' ...
  'movieDir' ...
  'movieKind' ...
  'movieName' ...
  'numphaseslices' ...
  'offByOne' ...                    % fill2seg param
  'outprefix' ...                   % old segmoviephase,manualcheckseg param
  'override' ...                    % old manualcheckseg param
  'partialDir' ...
  'partialfilename' ...             % old segmoviephase param
  'quickMode' ...                   % compileschnitz param
  'radius' ...
  'rback' ...                       % compileschnitz param
  'regsize' ...                     % old segmoviephase,manualcheckseg param
  'regsize' ...
  'rootDir' ...
  'schnitzName' ...                 % compileschnitz param
  'schnitzNum' ...                  % schnitzedit param
  'segall' ...                      % old segmoviephase param
  'segRange' ...
  'segmentationDir' ...
  'segoverride' ...                 % old segmoviephase param
  'segskip' ...                     % old segmoviephase param
  'trackMethod' ...                 % trackcomplete param
  'trackName' ...                   % trackcomplete param
  'trackRange' ...                  % trackcomplete param
  'trackUnCheckedFrames' ...        % trackcomplete param
  'transMax' ...                    % trackcomplete param
  'tracksDir' ...
  'upend' ...                       % old manualcheckseg param
  'usepartialfile' ...              % old segmoviephase param
  'yback0' ...                      % compileschnitz param
  'yfpOffset' ...                   % compileschnitz param
  };

upfields = upper(fields);

upperGivenField = upper(fieldNameString);

for i = 1:length(fields)
  if (length(upperGivenField) == length(char(upfields(i))))
    if (upperGivenField == char(upfields(i)))
      standardFieldString = char(fields(i));
      return
    end
  end
end

disp(['Warning: schnitzfield did not recognize field named ',fieldNameString,]);
disp(['         schnitzfield allowing ',fieldNameString,' to pass unchnaged.']);

standardFieldString = fieldNameString;
return
