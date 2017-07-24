function standardFieldString = DJK_schnitzfield(fieldNameString)
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

%   'Dskip' ...                       % old manualcheckseg param
%   'angThresh' ...
%   'autoCFL' ...                     % compileschnitz param
%   'autoYFL' ...                     % compileschnitz param
%   'autoGFL' ...                     % compileschnitz param
%   'autoRFL' ...                     % compileschnitz param
%   'beadCFP' ...                     % compileschnitz param
%   'beadYFP' ...                     % compileschnitz param
%   'beadGFP' ...                     % compileschnitz param
%   'beadRFP' ...                     % compileschnitz param
%   'cback0' ...                      % compileschnitz param
%   'cell' ...                        % fillinator,fill2seg param
%   'dateDir' ...
%   'expandvalue' ...                 % old manualcheckseg param
%   'fillinatorFile' ...              % fillinator and fill2seg param
%   'fillinatorRange' ...             % fillinator param
%   'finetuneimage' ...               % old manualcheckseg param
%   'flatFieldFile' ...               % compileschnitz param
%   'frnum' ...                       % old manualcheckseg param
%   'gback' ...                       % compileschnitz param
%   'imNumber1' ...
%   'imNumber2' ...

fields = {...
  'edge_lapofgauss_sigma' ...
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
  'prettyPhaseSlice' ...            % DJK 071122
  'segmentationPhaseSlice' ...      % DJK 071122
  'DJK_saveDir' ...                 % DJK 071122
  'useBigMask' ...                  % DJK 071126
  'maxCellWidthConservative' ...    % DJK 071127
  'noKinky' ...                     % DJK 071127
  'onScreen' ...                    % DJK 071210
  'singleCellMode' ...              % DJK 080704
  'fitRange' ...                    % DJK 080922
  'cropLeftTop' ...                 % DJK 081023
  'cropRightBottom' ...             % DJK 081023
  'edgeSlices' ...                  % DJK 081228
  'fillingEdgeSlices' ...           % DJK 081228
  'botHatSize' ...                  % DJK 081229
  'xlim' ...                        % DJK 090205
  'ylim' ...                        % DJK 090205
  'fieldColor_func' ...             % DJK 090505
  'lengthField' ...                 % DJK 090505
  'autoFluor' ...                   % DJK 090505
  'randomize' ...                   % DJK 090507
  'phaseImage' ...                  % DJK 090508
  'useArrayData' ...                % DJK 090508
  'binSize' ...                     % DJK 090508
  'selectionName' ...               % DJK 090508
  'fitTime' ...                     % DJK 090511
  'fluorShift' ...                  % DJK 090511
  'deconv_func' ...                 % DJK 090511
  'cropName' ...                    % DJK 090513
  'fitLength' ...                   % DJK 090517
  'plotRegression' ...              % DJK 090517
  'dataFields' ...                  % DJK 090518
  'normFields' ...                  % DJK 090518
  'sameLength' ...                  % DJK 090519
  'fieldY' ...                      % DJK 090519
  'fieldColor' ...                  % DJK 090519
  'colorlim' ...                    % DJK 090519
  'weightField' ...                 % DJK 090524
  'removeSmallCells' ...            % DJK 090604
  'frameSizes' ...                  % DJK 090605
  'lengthFields' ...                % DJK 090605
  'stabilize' ...                   % DJK 090618
  'colorbar' ...                    % DJK 090618
  'colorNormalize' ...              % DJK 090618
  'markerSize' ...                  % DJK 090618
  'fluor1' ...                      % NW 111202
  'fluor2' ...                      % NW 111202
  'fluor3' ...                      % NW 111202
  'weighing' ...                    % NW 120511 DJK_plotcrosscorr_stderror_store
  'bias'  ...                       % NW 120511 DJK_plotcrosscorr_stderror_store
  'extraNorm'  ...                  % NW 120511 DJK_plotcrosscorr_stderror_store
  'timeField'     ...               % NW 120511 DJK_plotcrosscorr_stderror_store
  'spacingError'      ...           % NW 120511 DJK_plotcrosscorr_stderror_store
  'dumpPlot' ...                    % MW 140428
  'setup' ...
  'softwarePackage' ...
  'camera', ...
  'customColors', ...               % MW 2015-02  
  'fluorcolor',...                  % MW 2015-04
  'minimalMode',...
  'rescaleCorrection',...           % MW 2015-04
  'dontOverwrite',...               % MW 2015-06
  'slices',...                      % MW 2015-06
  'rangeFiltSize',...
  'maskMargin',...
  'LoG_Smoothing',...
  'GaussianFilter',...
  'minDepth',...
  'neckDepth',...
  'assistedCorrection',...
  'showAll',...
  'rzeroonly'...
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

disp(['Warning: DJK_schnitzfield did not recognize field named ',fieldNameString,]);
disp(['         DJK_schnitzfield allowing ',fieldNameString,' to pass unchanged.']);

standardFieldString = fieldNameString;
return
