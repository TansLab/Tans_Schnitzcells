function slowSchnitzes=NW_detectSlowSchnitzes(p,schnitzcells,muField,varargin)
% finds schnitzes that grow very slowly and returns their numbers. Can be
% used as suggestions to delete schnitzes from analysis.
% Uses muField (which must be a growth rate) to test for slow growth
% Extension 2018-04: also mark those that have an imaginary component as
% slow.
% 
% Output: -  array with schnitznumbers
%         -  writes also a file in /analysisDir/slowschnitzes where detailed
%            frames and growth rates are stored for each schnitz
%
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzcells'  - schnitzcells structure
% 'muField'       - growth rate. Typically muP[xx]_fitNew where [xx] is a
%                   number that varies
%                   Choose one with averaged mu's!!
% 
% OPTIONAL ARGUMENTS:
% muThreshold     - threshold below which schnitzes are considered
%                   superslow growing. Default: 0 (i.e. all detected
%                   Schnitzes did shrink in at least one frame)


%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 3; functionName = 'NW_detect_slow_schnitzes';

if (nargin < numRequiredArgs) | (mod(nargin,2) ~= (mod(numRequiredArgs,2)) | ~isSchnitzParamStruct(p))
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
    p.(fieldName) = varargin{i+1};
  end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Parse the input arguments
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% set threshold for slow-growth detection
if ~existfield(p,'muThreshold ')
  p.muThreshold=0;
  disp('Treshold set to default value (0).');
end
% Set saveDir
if ~existfield(p,'NW_saveDir')
  p.NW_saveDir = [p.analysisDir 'slowschnitzes' filesep];
end
% Make sure that NW_saveDir directory exists
if exist(p.NW_saveDir)~=7
  [status,msg,id] = mymkdir([p.NW_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Not plotting!');
  return;
end

if ~existfield(schnitzcells(1),muField)
  disp(['Field ' muField ' does not exist. Not plotting!']);
  return;
end
%--------------------------------------------------------------------------
% Get fluorescence frames (muP[xx]_fitNew is only defined for these)
%--------------------------------------------------------------------------
fluor_frames=genvarname([upper(p.fluor1) '_frames']);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Open file to write schnitz numbers, growth rates and frame numbers.
%--------------------------------------------------------------------------
fid = fopen([p.NW_saveDir  p.movieName '-slowschnitzes.txt'],'wt');
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Checking for slowly growing Schnitzes.']);
dispAndWrite(fid, ['Growth rate variable is: ' muField]);
dispAndWrite(fid, ['Growth rate threshold is ' num2str(p.muThreshold)]);
dispAndWrite(fid, ['-------------------------------------------------']);   


%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Check Schnitzcells
%--------------------------------------------------------------------------
% array with slowly growing schnitzes
slowSchnitzes=[];

for schnitzrun=1:length(schnitzcells)
    
      isBelowThreshold=find(schnitzcells(1,schnitzrun).(muField)<p.muThreshold);
      
      % Also detect NaN values
      isBelowThreshold=[isBelowThreshold, find(isnan(schnitzcells(1,schnitzrun).(muField)))]; % MW 2015/07
      % Also detect imaginary values
      isBelowThreshold=[isBelowThreshold, find(   imag(schnitzcells(1,schnitzrun).(muField)) > 0  )]; % MW 2018/04
      
      if ~isempty(isBelowThreshold); 
          slowSchnitzes=[slowSchnitzes;schnitzrun]; % add schnitz to output-array
          % keep track of frames and growth rate in output file
          dispAndWrite(fid, ['Schnitz ' num2str(schnitzrun) ':']);
          dispAndWrite(fid, ['frame         mu (' muField ')']);
          for i=isBelowThreshold
              eval(['dispAndWrite(fid, [num2str(schnitzcells(1,schnitzrun).' fluor_frames '(i),3)     ''           '' num2str(schnitzcells(1,schnitzrun).(muField)(i),3)]);'])
              
          end
          dispAndWrite(fid,['...............'])
      end
end

% close file
fclose(fid);

disp('Identified following as slow/NaN:');
disp(mat2str(slowSchnitzes));





