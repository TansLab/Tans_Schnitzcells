function NW_saveSchnitzcells(p,myschnitzcells,varargin)
% This function saves a schnitzcells structure. Useful after rate addition
% and cell cycle corretion because there schnitzcells is not automatically
% saved.
% ************************ !! ******************
% Per default, myschnitzcells is renamed and saved as variable name 'schnitzcells' in a
% .mat file with name 'posXcrop-Schnitz.mat' (thus overwriting the standard
% saved schnitzcells structure).
% **********************************************
%
% ---------- REQUIRED ARGUMENTS -----------
% p
% myschnitzcells
% 
% ---------- OPTIONAL ARGUMENTS -----------
% 'NW_saveDir'        save directory. Default is p.tracksDir
%                     (../../posXcrop/data/). (Use filesep at end!)
% 'schnitzcellsName'  name of schnitzcells structure. Default is
%                     'schnitzcells'
% 'matFileName'       name of .mat file in which schnitzcells are stored.
%                     Default is [p.movieName '-Schnitz.mat]
%                     (posXcrop-Schnitz.mat)


%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
%% Settings
numRequiredArgs = 2; functionName = 'DJK_addToSchnitzes_length';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
  errorMessage = sprintf('%s\n%s',['Error width input arguments of ' functionName],['Try "help ' functionName '".']);
  error(errorMessage);
end

if ~isstruct(myschnitzcells)
  error('Second argument must be a schnitzcells structure.');
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
  for i=1:2:(numExtraArgs-1)
    if (~isstr(varargin{i}))
      errorMessage = sprintf('%s\n%s',['This input argument should be a String: ' num2str(varargin{i})],['Try "help ' functionName '".']);
      error(errorMessage);
    end
    fieldName = DJK_schnitzfield(varargin{i});
    p.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% Set default parameter values if they don't exist yet
% Save directory
if ~existfield(p,'NW_saveDir')
  p.NW_saveDir = [p.tracksDir];
end
% Make sure that NW_saveDir directory exists
if exist(p.NW_saveDir)~=7
  [status,msg,id] = mymkdir([p.NW_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.NW_saveDir ' : ' msg]);
    return;
  end
end
if ~existfield(p,'schnitzcellsName')
  p.schnitzcellsName = 'schnitzcells';
end
if ~existfield(p,'matFileName')
  p.matFileName = [p.movieName,'-Schnitz.mat'];
end
% create complete path for saving of .mat file
p.schnitzName = [p.NW_saveDir p.matFileName];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Save schnitzcells - Todo: Include check, whether saving will overwrite a
% former version of schnitzcells
%--------------------------------------------------------------------------
% rename variable
command=[p.schnitzcellsName '=myschnitzcells;'];
eval(command); %e.g. schnitzcells=myschnitzcells
command = ['save(p.schnitzName, ''' p.schnitzcellsName ''');'];
eval(command); % e.g. save(p.schnitzName, 'schnitzcells')

disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------

