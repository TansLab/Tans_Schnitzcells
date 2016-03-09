% DJK_selSchitzesToPlot selects a set of schnitzes from schnitzcells, based
% on provided criteria. Selected Schnitzes get field useForPlot=1, others
% useForPlot=0. If useForPlot is already set in provided schnitzcells,
% schnitzes with useForPlot=0 cannot become 1, unless the optional argument
% multiSel is set to 'OR'.
%
% If field is 'schnitzNr', will use schnitzNr (is not a field in
% schnitzcells, but the index)
%
% OUTPUT
% 's'               schnitzcells structure, with useForPlot added
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'    schnitzcells structure, in which useForPlot will be set
% 'field'           schnitzcells field, that will be used for criteria
% 'func'            function that specifies criteria (eg. @(x) x < 1)
%
% OPTIONAL ARGUMENTS:
% 'manualRange'     Allows to analyze a subset of frames
% 'ANDselection'=0  Either do OR or AND on subsequent selections (default: 1 = AND)
%

function s = DJK_selSchitzesToPlot(schnitzcells, field, func, varargin);

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 3; functionName = 'DJK_selSchitzesToPlot';
in = struct;

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
    in.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
if ~existfield(in,'manualRange')
  manualRange = false; % will use all cells
else
  manualRange = true;
end

if ~existfield(in,'ANDselection')
  in.ANDselection = 1; 
end
ANDselection = 1; % AND setting
if (in.ANDselection==0), ANDselection = 0; end % OR setting
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% In case useForPlot has not been set, set all to 1 (AND) or 0 (OR)
%--------------------------------------------------------------------------
s = schnitzcells;
if ~existfield(s(1), 'useForPlot')
  for cell = 1:length(s)
    % in case AND initially plotted
    % in case OR  initially not plotted
    s(cell).useForPlot = ANDselection; 
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% OR OPTION
%--------------------------------------------------------------------------
for cell = 1:length(s)
  % first check whether cell is in p.manualRange, if not continue
  if manualRange 
    cellFrames = s(cell).frames-1;
    if length(intersect(extraInputs.manualRange, cellFrames)) ~= length(cellFrames)
      if (ANDselection), s(cell).useForPlot = 0; end
      continue;
    end
  end
  
  if strcmp(field, 'schnitzNr') % use schnitzNr as field
    if func( cell )
      if (~ANDselection) s(cell).useForPlot = 1; end
    else
      if (ANDselection), s(cell).useForPlot = 0; end
    end
  else
    if func( s(cell).(field) )
      if (~ANDselection) s(cell).useForPlot = 1; end
    else
      if (ANDselection), s(cell).useForPlot = 0; end
    end
  end
end
%--------------------------------------------------------------------------



% % completeCycle
% s = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) x =~ 0, 'D', @(x) x =~ 0, 'E', @(x) x =~ 0);
% 
% % cells with data for more than 30 minutes
% s = DJK_selSchitzesToPlot(schnitzcells, 'time', @(x) x(end)-x(1) > 30);
% 
% % cells that are more than 1 um away from colony edge
% s = DJK_selSchitzesToPlot(schnitzcells, 'distCellToEdge', @(x) min(x)>1);
