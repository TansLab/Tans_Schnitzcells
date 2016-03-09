function schnitzcells = NW_addToSchnitzcells_noise(p, schnitzcells, varargin);
% NW_addToSchnitzcellsNoise adds the noise to each given dataField.
%
% Definition noise: normalized by subtracting mean for this timepoint off all branches.
% (norm term is not added since never used)
%
% Usually, noise is added with the function DJK_addToBranches_noise after branches 
% were extracted with DJK_getBranches. This is more direct since the i'th
% entry of each branch corresponds to time point i. 
% However, for NW_plot_cellCycleDependence data is needed in 'schnitzcells'
% and not 'branches' structure, since a time trace for each single schnitz will 
% be drawn.
% NW_addToSchnitzcells_noise internally calls DJK_getBranches and
% DJK_addToBranches_noise and then writes the calculated noise for each
% schnitz and time point back into the schnitzcells structure.
%
% There is probably a more direct way to calculate the mean of the
% dataField at each time point (if, and  I think so, the mean is calculated
% via the arithmetic mean without extra weighing of branches), however,
% the used way will avoid errors and make sure, that the added noise terms
% are identical to the ones used in correlation plots. (NW 2012-03)
%
%
% OUTPUT
% 'schnitzcells'     - schnitzcells structure with added noise terms
%                      SO FAR THIS SCHNITZCELLS_NOISE IS NOT SAVED but only
%                      present in the workspace
%
% REQUIRED ARGUMENTS:
% 'p'
% 'schnitzcells'      schnitzcells structure
%
% OPTIONAL ARGUMENTS:
% 'dataFields'        fields to be stored in branches
%                     default: {'Y_time' 'Y6_mean' 'muP15_fitNew'}
%                     I think the first entry must be a Time. Always choose
%                     a time that appears as frequently as the variable!!
%                     ('time' vs 'Y_time')
% 'noisyBranchData'   branch data where noise has already been added 
%                     (don't restrict time interval and don't use only cycle cells! 
%                     take all schnitzes into account that you want to plot later!). If 
%                     not given, branchdata+noise will be calculated directly here.
%                     **** 
%                     Should only be used for quick debugging since risk of
%                     giving not corresponding schnitzcells+branchdata
%                     *****
%
% if you want to add noise to data, that appears in a different frequency
% (e.g. muP15_fitNew_all vs Y6_mean) you have to run the program twice with
% respetive time fields
%
% *****************************************************************
% Optimally, parse schnitzcells, that are neither restricted to fitTime nor
% cyclic cells (e.g. 's_rm'). If not (and maybe also if yes), many branches
% will have different lengths and by default schnitzes that correspond to
% shorter branches are ignored when adding noise.
% You can avoid this problem by changing 'sameLength' option in
% DJK_getBranches below. (I am however unsure if then noise terms are
% slightly different) (NW 2012-03)%
% *****************************************************************
%

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 2; functionName = 'NW_addToSchnitzcells_noise';

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
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If not provided, use standard dataFields
if ~existfield(p, 'dataFields')
  p.dataFields = {'Y_time' 'Y6_mean' 'muP15_fitNew'};
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Check if dataFields is empty. Then no noise has to be added
%--------------------------------------------------------------------------
% If not provided, use standard dataFields
if isempty(p.dataFields)
    disp(['No fields chosen to which noise should be added. Exiting...']);
    return
end
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% Check Schnitzcells
%--------------------------------------------------------------------------
if length(schnitzcells) < 1
  disp('Schnitzcells is empty. Not plotting!');
  return;
end

%if dataFields have different length/frequency -> cannot work
for i = 2:length(p.dataFields)
    if length(schnitzcells(1).(char(p.dataFields(2))))~=length(schnitzcells(1).(char(p.dataFields(1))))
        disp(['dataFields have different frequency. Exiting...']);
            return
    end
end
% Existence of DataFields will be checked directly in DJK_addToBranches_noise
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get BranchData and add Noise to BranchData. Create Noise names
%--------------------------------------------------------------------------
if ~existfield(p, 'noisyBranchData')
    schnitzcells;
    myBranchData = DJK_getBranches(p,schnitzcells,'dataFields',p.dataFields,'sameLength',0); % don't restrict fitTime!
    %myBranchData = DJK_getBranches(p,schnitzcells,'dataFields',p.dataFields); % don't restrict fitTime!
    myNoisyBranchData = DJK_addToBranches_noise(p, myBranchData,'dataFields',p.dataFields);
else
    myNoisyBranchData = p.noisyBranchData;
end
    
for i = 2:length(p.dataFields)
  field = char(p.dataFields(i));
  noisefield{i} = ['noise_' field];
  %normfield{i}  = ['norm_' field]; not used
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ASSOCIATE NOISE TO EACH SCHNITZ IN SCHNITZCELLS
%--------------------------------------------------------------------------
% if a schnitz appears in several branches, it does not matter from which
% branch the noise data is taken since it is always the same.
%              -----schn2----- (branch1)
%            /
% -schn1-----
%            \
%              -----schn3----- (branch2)
% noise terms of schn1 are the same in branch1 and branch2


%--------------------------------------------------------------------------
% loop over all schnitzes
%--------------------------------------------------------------------------
for schnitz = 1:length(schnitzcells)
    % let know how far we are
    if mod(schnitz,100)==0
        disp(['Adding noise to schnitz ' num2str(schnitz) '...']);
    end
    %search schnitz in all branches
    for branchrun=1:length(myNoisyBranchData)
        schnitz_idx=find(myNoisyBranchData(branchrun).schnitzNrs==schnitz);
        %if schnitz exists in branch
        if ~isempty(schnitz_idx)
            %check if #frames and exact number of frames in which schnitz
            %appears in branch data is identical to schnitzcell info
            %(should be)
            
            %debugging
            %schnitz=schnitz;
            %branchrun=branchrun;
            %schnitz_idx=schnitz_idx;
            %myNoisyBranchData(branchrun).(char(p.dataFields(1)));
            %schnitzcells(schnitz).(char(p.dataFields(1)));
            
            branchTimes=myNoisyBranchData(branchrun).(char(p.dataFields(1)))(schnitz_idx); %should be identical
            schnitzcellTimes=schnitzcells(schnitz).(char(p.dataFields(1)));
            
            % number of frames wrong 
            if length(schnitz_idx)~=length(schnitzcells(schnitz).(char(p.dataFields(1))))
                disp(['Number of frames in which schnitz ', num2str(schnitz), ' appears is inconsistent. Exiting...']);
                return
                
            % not all frame numbers are identical (i.e. different time
            % points)
            elseif (sum(branchTimes==schnitzcellTimes)~=length(schnitz_idx) ) 
                disp(['Specific frame numbers of schnitz ', num2str(schnitz), ' are inconsistent. Exiting...']);
                return
                
            % everything allright. Add noise info, then continue with next
            % schnitz
            else
                % loop over all dataFields (except time)
                for i=2:length(p.dataFields)
                    schnitzcells(schnitz).(noisefield{i})=myNoisyBranchData(branchrun).(noisefield{i})(schnitz_idx);
                end
                continue
            end
        end
    end
end
        



