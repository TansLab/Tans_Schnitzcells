function p = PN_manualcheckseg (p, varargin);
% Lc will be saved as the corrected image
%
%   MANUALCHECKSEG allows users to review the results of image segmentation
%
%   in order to manually correct any merged cells or delete false positives.
%
%   MANUALCHECKSEG(P,'Field1',Value1,'Field2',Value2,...) also performs
%   manual/interactive correction of segmentation results, but permits users
%   to adjust any parameters describing the movie or parameters controlling
%   the manual segmentation checking process.  The routine first sets all of
%   the movie analysis parameters to their default values defined for the given
%   movieKind, and then overwrites any specific parameters provided by setting
%   P.Field1 = Value1, P.Field2 = Value2, etc.  Thus any/all schnitzcells
%   parameter values can be defined in the function call via these optional
%   field/value pairs.  (This is in the style of setting MATLAB properties
%   using optional property/value pairs.)
%
%   MANUALCHECKSEG returns a struct (1x1 struct array) referred to as the
%   schnitzcells parameter structure that contains fields and values
%   contained in P, including unchanged/original parameters plus any of those
%   added or overridden in the list of properties & values provided in the
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control manualcheckseg
%
%   outprefix     overwrites [p.movieName 'seg']
%   manualRange   specify frame range to check
%   overwrite     if 1, program will redo frames already having Lc; default=0
%   expandvalue   increases each segmentation image border by this amount
%                 (in case the segmentation cut off some cells); defaults to 30
%   frnum         check one frame only, also sets p.overwrite = 1
%   Dskip         determines how many frames are skipped with e.g. 
%                 '.' and ',' commands.
%   upend         frame size / location; default is 680
%   leftend       frame size / location; default is 516
%   min_size      [height width] of figure; default [upend-60 leftend-10]
%   finetuneimage if 1, program will not renumber the image; default = 0
%   regsize       maximum size in pixels of any translation between phase and
%                 fluorescent images; default is 3
%   assistedCorrection      displays potentially wrong segmented cells in
%                           (added by PN the 23-04-2012)
%                    Works sometimes only if every frame is used (or stepsize not too big)! otherwise error
%                    when backing up. exact conditions...no idea (comment NW2012-05-10)
%   maxImage      if 1, image will be resized to fill screen. default=1
%   showAll       default =1; if showAll=1 will show all frames to user.
%                 (not sure of its behavior, testing still -MW 2015/06)
%   problemCells  if 'problemCells' (or p.problemCells) is set, then cells
%                 that are given by this array, which should consist of 
%                 [schnitz, framenr, cellnr; etc..] problem cells will be
%                 highlighted too.
%
%-------------------------------------------------------------------------------
%

%-------------------------------------------------------------------------------
% variable renaming during rewrite:
%   outprefix -> p.outprefix
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------


% elapsed time for frame 444 in 2012-05-08. 390 cells


numRequiredArgs = 1;
if (nargin < 1) | ...
        (mod(nargin,2) == 0) | ...
        (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> manualcheckseg:',...
        '    Invalid input arguments.',...
        '    Try "help manualcheckseg".');
    error(errorMessage);
end

if ~existfield(p,'showAll')  
    p.showAll=1;  
end



% Load the watermark (MW edit 2014/12)
% Note that only the green channel is used, and only binary form.
p.mywatermark=imread('watermark.png');

%-------------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for frameIdx=1:2:(numExtraArgs-1)
        if (~isstr(varargin{frameIdx}))
            errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
                'Error using ==> manualcheckseg:',...
                '    Invalid property ', num2str(varargin{frameIdx}), ...
                ' is not (needs to be) a string.',...
                '    Try "help manualcheckseg".');
            error(errorMessage);
        end
        fieldName = DJK_schnitzfield(varargin{frameIdx}); % DJK 090210
        p.(fieldName) = varargin{frameIdx+1};
    end
end


disp(' ')
disp('Instructions: press left mouse button (or key ''4'') on two consecutive cells to join them.')
disp('              Press right mouse button (or key ''5'') to cut at that point.')
disp('              Shift + left mouse button (or key ''7'') erases current cell.')
disp('              When you''re done, press <space> to go to the next frame.')
disp('              press <escape> to redo this frame from the original autoseg file.')
disp('              to save progress during work on a frame, press:')
disp('                    ''w'' to write a temporary partial correction to the file.')
disp('              press ''x'' to black out an area.')
disp('              press ''t'' to mark terraced area.')
disp('              press ''v'' to define region of interest and restrict to this selection.')
disp('              press ''c'' to morphologically close each cell area (imclose)')
disp('              press ''o'' to obliterate all but the cell you''re pointing to.')
disp('              press ''.'' to skip frame (without saving).')
disp('              press '','' to go back a frame (without saving).')
disp('              press ''a'' to add a new cell where the mouse is pointing.')
disp('              press ''q'' to quit.')
disp(' ')
disp('              press ''e'' to expand image.')
disp('              press ''f'' for "fine-tuning" (to avoid renumbering the image).')
disp('              press ''g'' to goto indexnum = ... .')
disp('              press ''r'' to show 1µm long bar.')
disp('              press ''i'' to fill the cell you are pointing to.')
disp('              press ''j'' to fill all cells.')
disp('              press ''h'' to reseed a cell.')
disp('              press ''d'' to remove all ''cells'' which look like dirt (large convex hull).')
disp('              press ''u'' to go back one step (possible once).')
disp('              press ''m'' to switch between phase contrast and segmentation image.')
disp('              press ''l'' to show perimeters of cell areas.')
disp('              press ''b'' to restrict segmentation to overlaps with previous segmentation image')
disp('                             (requires assistedCorrection mode, protential trouble after backstep '','')')
disp('              press ''n'' to toggle fullscreen mode (go to next frame after toggling it off).')
disp('                    ''s'' to toggle no numbers, cell numbers, schnitz numbers.')
disp('              press ''k'' to merge cells in previous frame that have joined into one cell in this frame')

disp(' ')

%%
quit_now=0;
global pos Limage ourfig pp phfig showPhase figureToFocusOn % flfig % DJK 071206

% create new figures
%phfig  = figure(); clf; % MW should be done later when needed
ourfig = figure(); clf;

% phfig/ourfig settings
figureToFocusOn = ourfig;
showPhase=0;

% flfig  = figure; % former version: fluor picture
outl=length(p.segmentationDir);

% set some defaults if they don't exist
if ~existfield(p,'Dskip')
    p.Dskip=1;
end

if ~existfield(p,'outprefix')
    p.outprefix = [p.movieName 'seg'];
end

D = dir([p.segmentationDir, p.outprefix, '*.mat']);

if isempty(D)
   error('Could not find files, mayb be path not set correctly?'); 
end

[s,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;
if ~existfield(p,'manualRange')
    % p.manualRange=1:length(dir([p.segmentationDir, p.outprefix, '*.mat']));
    segNameStrings = char(s);
    p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
end

if ~existfield(p,'overwrite')
    p.overwrite=0;
end
if existfield(p,'override') % backward compatibility; legacy
    disp('Please use p.overwrite instead of p.override');
    p.overwrite=p.override;
end

if ~existfield(p,'expandvalue'),
    p.expandvalue = 30;
end;

if existfield(p,'frnum')
    %  p.manualRange=1:length(dir([p.segmentationDir, p.outprefix, '*.mat']));
    newloopindex = find(p.manualRange==frnum);
    if isempty(newloopindex)
        disp(['Could not check frame ', num2str(gotoframe), ...
            ' because it is not in the manualRange of segmented files.']);
    else
        p.manualRange=frnum;
        p.overwrite=1;
    end
end

if ~existfield(p,'screensize')
    try
        p.screenSize=get(0,'ScreenSize');
    catch
        p.screenSize=[1280 1024];
        warning('Set screen size to manual dimensions');
    end
end

if ~existfield(p,'min_size') %[height width] desired figure, using 70 for top part, 35 for bottom part, 10 for sides
    p.min_size = p.screenSize(3:4)-[105 20]; % DJK 090117 was [p.upend-60 p.leftend-10];
end

if ~existfield(p,'finetuneimage')
    p.finetuneimage = 0;
end

if ~existfield(p,'figs')
    p.figs = 0;
end

if ~existfield(p,'fill_cut')
    p.fill_cut = 1;
end

if ~existfield(p,'regsize')
    p.regsize = 3;
end

if ~existfield(p,'assistedCorrection')
    p.assistedCorrection = 0;
end
if p.assistedCorrection
    disp('You have chosen the Assisted Correction Mode :')
    disp(' - potentially re-merged or too small cells are displayed in white;')
    disp(' - centroids of cells of the former image are indicated as black dots.')
end
if ~existfield(p,'maxImage')
    p.maxImage=1;
end

% get screenSize for image display
myScreenSize=get(0,'ScreenSize');
maxValidImageSize=myScreenSize(3:4)-[150,150]; % empirical
%


backwards = 0;
gotoframenum=0;
loopindex = 1;

%% main loop
while loopindex <= length(p.manualRange)
   frameIdx = p.manualRange(loopindex);
   
   %former image data if available (L_prec for L preceeding)
    L_prec=[];
    rect_prec=[];
    if frameIdx > p.manualRange(1) %if possible loads the preceeding segmented image, this is used only for assisted correction
        
        filename = [p.segmentationDir,p.movieName,'seg',str3(p.manualRange(loopindex-1)),'.mat']; % MW 2015/08
        
        if exist(filename,'file')
            load(filename);
            if exist('Lc','var')
                 L_prec=Lc;
                 rect_prec=rect;
            end
            clear Lc LNsub rect
        end
        
    end
    
    %new image data
    clear Lc creg yreg savelist rect newrect oldrect;
    name= [p.segmentationDir,p.movieName,'seg',str3(frameIdx)];
    tempsegcorrect=0;
    load(name);
    
    % nitzan's changes June24th
    if exist('phaseFullSize')==1
        p.fullsize = phaseFullSize;
    end
    
    if ~isfield(p,'fullsize')
        p.fullsize = [1024 , 1344];
        disp('setting size of images arbitrarily to [1024 , 1344].');
    end
    
    %{
    if (gotoframenum && exist('Lc')==1 && p.overwrite==false)
        p.overwrite=true; % MW todo ?! 
        disp('WARNING: p.overwrite=true (MW todo; doesn''t seem correct behavior).');
    end %DJK 071207
    %}
    
    % Set flag if file was already checked before.
    if exist('Lc')==1
        LNsub=Lc;
        p.CurrentFrameApprovedFlag = 1; % Addition MW 2014/12
    else
        p.CurrentFrameApprovedFlag = 0; % Addition MW 2014/12
    end    
    
    % Parse current frame index
    p.currentFrame = frameIdx;
    
    % nitzan's changes June24th
    if p.showAll || ... % MW 2015/06
        ((exist('Lc')~=1 || backwards || p.overwrite==1 || (exist('Lc')==1 && tempsegcorrect==1)) && ~mod(frameIdx-1,p.Dskip))
        
        %% Calculate phase image
        g = double(phsub);
        if length(g) > 0
            g = DJK_scaleRange(g, [max(max(g)) min(min(g))], [0 1]);
            g = DJK_scaleRange(g, [0.2 1], [0 1]);
        end
        g = uint8(g * 255);            
        
        is_done=0;
        savelist=['''Lc'''];
        crop_pop=0;

        % DJK: multiple settings in 1 variable
        clear DJK_settings
        DJK_settings.finetuneimage = p.finetuneimage; DJK_settings.figs = p.figs; DJK_settings.fill_cut = p.fill_cut;
        
        %set(0,'CurrentFigure',ourfig)
        while ~is_done
           
            % This executes the function that asks for keyboard input.
            [p,Lc,is_done,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,gotoframenum,DJK_settings] = ...
                PN_manual_kant(p, LNsub, L_prec, g, rect,rect_prec,phsub, DJK_settings,p.assistedCorrection); % p added 2014/12 MW
                        
            if backwards && (loopindex>1),

                loopindex = loopindex-p.Dskip;

                disp(['backing up to frame = ',...
                    str3(p.manualRange(loopindex)),...
                    '(loopindex = ', num2str(loopindex), ')']);
                % JCR: Note this code is still kinda bad because when loopindex = 1
                % and user presses backwards, loopindex skips forward anyway.
                % Probably should simplify and get rid of Dskip all together,
                % now that segRange and manualRange should handle arbitrary frames
            end;
            
            if quit_now                                
                
                close(ourfig);
                if showPhase close(phfig); end;
                clear global pos Limage ourfig res pp phfig showPhase;                
                
                disp('Will return now. Bye.');                
                return;
                
            end;
            
            
            if savetemp
                
                LNsub=Lc;
                if savetemp==2
                    tempsegcorrect=1;
                    eval(['save(''',name,''',',savelist,',''tempsegcorrect'',''-append'');']);
                    disp(['Saved partial file ',name(outl:end),'   # of cells: ',num2str(max2(Lc))]);
                end
                
            end
            
        end
        
        
        % DJK: restore settings from 1 variable
        p.finetuneimage = DJK_settings.finetuneimage; p.figs = DJK_settings.figs; p.fill_cut = DJK_settings.fill_cut;
        
        if ~dontsave,
            
            tempsegcorrect=0;
            eval(['save(''',name,''',',savelist,',''tempsegcorrect'',''-append'');']);
            disp(['Updated file ',name(outl:end),'   # of cells: ',num2str(max2(Lc))]);
            
        elseif addtolist
            
            if exist(badseglistname)==2
                load(badseglistname);
            end
            if exist('badseglist')~=1
                badseglist=[];
            end
            badseglist=[badseglist,frameIdx];
            save(badseglistname,'badseglist');
            disp(['Added frame ',name(outl:end),' to bad segmentation list.']);
            
        else
            disp(['Skipped file ',name(outl:end),' <--']);
        end;
        
    end
    
    if ~backwards
        loopindex = loopindex + p.Dskip;                                     %here the loop index is incremented
    end
    
    if gotoframenum ~= 0
        
        newloopindex = find(p.manualRange==gotoframenum);
        
        if isempty(newloopindex)
            disp(['Could not goto frame ', num2str(gotoframenum), ...
                ' because it is not in the manualRange.']);
        else
            loopindex = newloopindex;
        end
        
    end
    
    %disp(['Now set loopindex to ', num2str(loopindex)]);
    
    
end

%%

close(ourfig);
if showPhase, close(phfig); end;
clear global pos Limage ourfig res pp phfig showPhase;

disp('bye');

end
