function p = NW_initializeFluorData(p, varargin)
%
% DOES NOT ACCOUNT FOR EXTRA RESCALE CORRECTION (WRONG SCALING OF CFP
% IMAGES). Influences data in creg(!), cback(?). Data will be overwritten
% anyways in DJK_correctFluorImage_anycolor
% FURTHERMORE ASSUMES PHASE IMAGE BINNING = 1; binning of the
% fluor images is inferred from the size ratio between the phase image and 
% the fluor image. However, if binning was also used for the phase image
% (typically should not be the case), then the binning will be determined
% incorrectly, which has consequences for the fluor correction later on.
% Use the  p.resizePhase option in DJK_cropImages_3colors.m to address this
% issue.
%
%
% TODO: INCLUDE SAVEDIR PROPERLY OR LEAVE OUT! (now: automatically in
% segmantation folder written
%
% Adds information about the fluorescence Data to pos[x]cropseg[xxx].mat.
% Initially this was part of the segmentation step. Now, segmentation of
% phasecontrast pictures and fluorescence analysis are seperated.
% All fluorescence colors are added within this single function.
%
% -----------------------------------------------------------
% Required Input Arguments
% -----------------------------------------------------------
% p                      struct parameter
%
%
% -----------------------------------------------------------
% Optional Input Arguments
% -----------------------------------------------------------
%
% manualRange           if given, only these images will be segmented.
%                       Otherwise whole range of images will be analyzed.
%
% NW_saveDir            directory in which .mat files will be stored.
%                       default: [p.analysisDir 'segmentation\']. This
%                       overwrites the existing .mat files with new files
%                       that now include fluorescence data. Default option
%                       seems to be the most useful one.



%-------------------------------------------------------------------------------
%% Parse the input arguments, input error checking. Use inputParser in
% recent matlab releases.
%-------------------------------------------------------------------------------
numRequiredArgs = 1;
if (nargin < 1) | ...
        (mod(nargin,2) == 0) | ...
        (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> NW_initializeFluorData:',...
        '    Invalid input arguments.',...
        '    Try "help NW_initializeFluorData".');
    error(errorMessage);
end

numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
                'Error using ==> NW_initializeFluorData:',...
                '    Invalid property ', num2str(varargin{i}), ...
                ' is not (needs to be) a string.',...
                '    Try "help NW_initializeFluorData".');
            error(errorMessage);
        end
        fieldName = DJK_schnitzfield(varargin{i}); % DJK 071122
        p.(fieldName) = varargin{i+1};
    end
end


%--------------------------------------------------------------------------
% Checking or creation of directories
%--------------------------------------------------------------------------
% make sure every directory field has a trailing filesep
%if isfield(p,'NW_saveDir') & (p.NW_saveDir(end) ~= filesep)
%    p.NW_saveDir = [p.saveDir filesep];
%end
% if directory doesn't exist, create it
%if exist(p.NW_saveDir)~=7
%    [status,msg,id] = mkdir([p.segmentationDir p.NW_saveDir]);
%    if status == 0
%        disp(['Warning: unable to mkdir ' p.NW_saveDir ' : ' msg]);
%    end
%end


%-------------------------------------------------------------------------------
% Figure out which seg-files of images exist
%--------------------------------------------------------------------------

%  % figure out the number of phase image layers (slices) per position capture
%Dphase1   = dir([p.imageDir, '*-p-1-*.tif']);
%DphaseAll = dir([p.imageDir, '*-p-*.tif']);
%if ~isfield(p,'numphaseslices')
%    if isempty(DphaseAll)
%        error(['Can''t find any images in directory ' p.imageDir]);
%    else
%        p.numphaseslices = int8(ceil(length(DphaseAll)/length(Dphase1)));
%        disp(['You appear to have ' num2str(p.numphaseslices) ' phase image per frame.']);
%    end
%end

%if ~isfield(p,'slices')
%    p.slices = [ 1:p.numphaseslices ]; % by default, use all slices
%end

%% check if segmentedImages exist and are in folder /segmentation/
segmentedImages=dir([p.segmentationDir, '*seg*.mat']); % MW 2015/04 
if isempty(segmentedImages)
    error(['Can''t find any matlab files of segmented images in directory ' p.segmentationDir])
end


% if user didn't specify a range of frames for which to add fluor data,
% figure it out by counting ".mat" files of segmented frames in
% /segmentation/
if ~isfield(p,'manualRange')
    error('manualRange not given in optional parameters.');
        % NOTE: the code below would figure out the range by itself, but
        % does not work for numbers >999. This could be fixed by using more
        % advanced regular expression.
        % Therefor, let's just throw an error here now.
    [s,I] = sort({segmentedImages.name}');
    segmentedImages = segmentedImages(I);
    numpos= findstr(segmentedImages(1).name, '.mat')-3;
    imageNameStrings = char(s);
    p.manualRange = str2num(imageNameStrings(:,numpos:numpos+2))';
    clear s I numpos imageNameStrings;
end


%--------------------------------------------------------------------------
% Display settings that will be used
disp ('using schnitzcells parameter structure:');
disp (p);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Loop over frames : load .mat files of segmented Images and add fluor data
%--------------------------------------------------------------------------
% set fluorescence picture counters to '-1' (i.e. fluor color not
% attributed)
fluor1counter=-1; fluor2counter=-1; fluor3counter=-1;

for i= p.manualRange
    %%
    %----------------------------------------------------------------------
    % Search for [segmentation].mat file
    segImage =  dir([p.segmentationDir, '*seg', str3(i), '.mat']); % .mat file of current image % MW 2015/04
   
    % check if segmentation of this frame exists (and not falsely imposed
    % by manualRange)
    if isempty(segImage)
        message=sprintf('Cannot find segmented image with frame number %d. Will skip number...',i);
        disp(message);
    else
        %% Load phaseFullSize, LNsub and rect from [segmentation].mat file
        eval(['load(''' p.segmentationDir segImage.name ''', ''phaseFullSize'',''LNsub'',''rect'')']);
         
        %% ------------------extract data ---------------------------------
        %
        % NOTE that this code "loops" over the 3 fluor codes in a
        % hard-coded way, i.e. code has been copied-pasted three times, for
        % each possible fluor once. 
        % MW 2015/04
        %
        % fluor1: check for existence of image and add data to .mat file if
        % existent
        if ~strcmp(p.fluor1,'none')
            
            if (fluor1counter==-1) fluor1counter=0; end
            
            % filename of fluor data picture
            fluor1name=[p.imageDir,p.movieName,'-' p.fluor1 '-',str3(i),'.tif'];
            
            if exist(fluor1name,'file')==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error
                
                %%
                disp('found fluor1 image');
                
                fluor1counter=fluor1counter+1;
                
                % obtain fluor data
                [fluor1reg, fluor1shift, fluor1back, fluor1binning] = quicknoreg(LNsub,fluor1name,rect,0,phaseFullSize); 
                % obtain info (EXPosure time, importantly)
                [exptfluor1str, gainfluor1, exptfluor1] = imsettings(fluor1name);
                
                % change names to actual color: fluor1back -> yback, or
                % rback etc.. fluor1 could in principle be kept in the
                % variables. However, when e.g. falsely initiating later
                % with twisted fluor1 and fluor2 colors, the information
                % which real color the fluorescence picture corresponds to,
                % is lost.
                % not too beautiful but doesn't cost a lost of time ...
                reg=genvarname([p.fluor1 'reg']); eval([reg '=fluor1reg;']);
                %shift=genvarname([p.fluor1 'shift']); eval([reg '=fluor1shift']);
                back=genvarname([p.fluor1 'back']); eval([back '=fluor1back;']);
                binning=genvarname([p.fluor1 'binning']); eval([binning '=fluor1binning;']);
                expt=genvarname([ 'expt' p.fluor1]); eval([expt '=exptfluor1;']);
                gain=genvarname(['gain' p.fluor1]); eval([gain '=gainfluor1;']);
                % exptfluor1str ...
                
                % save to existing [segmentation].mat file (append)
                eval(['save(''' p.segmentationDir segImage.name ''',''-append'',''' reg ''',''' back ''',''' binning ''',''' expt ''',''' gain ''');']);
                clear reg back binning expt gain;
            else
                disp(['Info: ' fluor1name ', has no fluor..']);
            end
        end
                
        % fluor2: check for existence of image and add data to .mat file if
        % existent
        if ~strcmp(p.fluor2,'none')
            if (fluor2counter==-1) fluor2counter=0;end
            fluor2name=[p.imageDir,p.movieName,'-' p.fluor2 '-',str3(i),'.tif'];
            if exist(fluor2name)==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error
                disp('found fluor2 image');
                fluor2counter=fluor2counter+1;
                [fluor2reg, fluor2shift, fluor2back, fluor2binning] = quicknoreg(LNsub,fluor2name,rect,0,phaseFullSize); 
                [exptfluor2str, gainfluor2, exptfluor2] = imsettings(fluor2name);                                               
                
                % change names to actual color: fluor2back -> yback, or
                % rback etc.. 
                reg=genvarname([p.fluor2 'reg']); eval([reg '=fluor2reg;']);
                %shift=genvarname([p.fluor2 'shift']); eval([reg '=fluor2shift']);
                back=genvarname([p.fluor2 'back']); eval([back '=fluor2back;']);
                binning=genvarname([p.fluor2 'binning']); eval([binning '=fluor2binning;']);
                expt=genvarname([ 'expt' p.fluor2]); eval([expt '=exptfluor2;']);
                gain=genvarname(['gain' p.fluor2]); eval([gain '=gainfluor2;']);
                % exptfluor2str ...
                
                % save to existing [segmentation].mat file (append)
                command = ['save(''' p.segmentationDir segImage.name ''',''-append'',''' reg ''',''' back ''',''' binning ''',''' expt ''',''' gain ''');'];
                eval(command);
                clear reg back binning expt gain;
            end
        end
        
        % fluor3: check for existence of image and add data to .mat file if
        % existent
        if ~strcmp(p.fluor3,'none')
            if (fluor3counter==-1) fluor3counter=0;end
            fluor3name=[p.imageDir,p.movieName,'-' p.fluor3 '-',str3(i),'.tif'];
            if exist(fluor3name)==2 & numel(rect)>0 % if no cells found, rect is empty, and gives error
                disp('found fluor3 image');
                fluor3counter=fluor3counter+1;
                [fluor3reg, fluor3shift, fluor3back, fluor3binning] = quicknoreg(LNsub,fluor3name,rect,0,phaseFullSize); 
                [exptfluor3str, gainfluor3, exptfluor3] = imsettings(fluor3name);
                
                % change names to actual color: fluor3back -> yback, or
                % rback etc.. 
                reg=genvarname([p.fluor3 'reg']); eval([reg '=fluor3reg;']);
                %shift=genvarname([p.fluor3 'shift']); eval([reg '=fluor3shift']);
                back=genvarname([p.fluor3 'back']); eval([back '=fluor3back;']);
                binning=genvarname([p.fluor3 'binning']); eval([binning '=fluor3binning;']);
                expt=genvarname([ 'expt' p.fluor3]); eval([expt '=exptfluor3;']);
                gain=genvarname(['gain' p.fluor3]); eval([gain '=gainfluor3;']);
                % exptfluor3str ...
                
                % save to existing [segmentation].mat file (append)
                eval(['save(''' p.segmentationDir segImage.name ''',''-append'',''' reg ''',''' back ''',''' binning ''',''' expt ''',''' gain ''');']);
                clear reg back binning expt gain;
            end
        end
        % --------------------------------------------------------------
    end
end
% warn if maybe fluorescence color is set wrong
if (fluor1counter==0)
    disp('Found no Fluor1 image. Maybe you set the color wrong. Change with p.fluor1=...');
end
if (fluor2counter==0)
    disp('Found no Fluor2 image. Maybe you set the color wrong. Change with p.fluor2=...');
end
if (fluor3counter==0)
    disp('Found no Fluor3 image. Maybe you set the color wrong. Change with p.fluor3=...');
end
