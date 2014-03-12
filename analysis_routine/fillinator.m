function p = fillinator(p,varargin); % A new hope!

% FILLINATOR  Manually flood fill and track an individual cell across frames.
%
%   FILLINATOR allows users to track one cell identified by manual clicking
%   across a series of frames within a movie.  The cell in one frame is
%   implicitly matched to the selected cell in adjacent frames.  In the words
%   of the author (Jordi):
%
%     "... a small code that lets [users] go backwards in time,
%     tracking manually the cells and measuring the average YFP and
%     CFP. Following Gurol's suggestion, the code uses mainly the
%     fluorescence images for cell identification (which is done with
%     Michael's proposal: edge detection + dilation), even though phase
%     images can also be used. It takes me between 30min and an hour to
%     track a single competent event, and the measurement is much more
%     reliable than with the quick-and-dirty measure."
%
%   Look for interactive instructions output to the MATLAB command window.
%
%   FILLINATOR outputs a file known as the fillinatorFile that contains a 
%   struct array containing frame, xpixels, ypixels for each frame processed. 
%
%   FILLINATOR(P,'Field1',Value1,'Field2',Value2,...) does the same thing,
%   but the extra arguments permit users to adjust any parameters controlling
%   the process.  The extra arguments can override any specific parameters
%   provided in P by setting P.Field1 = Value1, P.Field2 = Value2, etc.
%   Thus any/all schnitzcells parameter values can be defined in the
%   function call via these optional field/value pairs.  (This is in the
%   style of setting MATLAB properties using optional property/value pairs.)
%
%   FILLINATOR returns a struct (1x1 struct array) referred to as the
%   schnitzcells parameter structure that contains fields and values
%   contained in P, including unchanged/original parameters plus any of those
%   added or overridden in the list of properties & values provided in the
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control FILLINATOR
%
%   cell             an integer indicating the number (ID) of a cell that is 
%                    used when naming the fillinator file
%   
%   fillinatorFile   filename of fillinator pixel data file generated (default
%                    filename is [p.segmentationDir 'cell-001-pixels.mat'])
%
%   fillinatorRange  a range of frame numbers on which to apply the fillinator
%
%   Note: Optional args e.g. 'cell',1 passed to initschnitz will define 
%         'p.cell' and 'p.fillinatorFile' for you and will append '-cell-001' 
%         to p.segmentationDir, p.tracksDir and p.partialDir, so that you can 
%         manage more than one fillinator cell per movie, each with it's own 
%         set of sub-directories.
%   
%   Note: In spite of specifying cell and fillinatorFile in P or in optional 
%         arguments, this function will still prompt the user to verify them.
%-------------------------------------------------------------------------------
%
% Actual operating instructions:
%    Once you're running, you can do the following:
%    Click on a cell to select it.
%    If you like it, hit return
%    Otherwise, press C, Y, or P to switch to cyan,yellow, or phase images
%    That's it.

% Programming Notes:
%
% 5/20/05-- i removed the dilat business, since we're not using it --ME
% 6/01/05-- added more cell,fillinatorFile code here and in initschnitz --JCR


%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
        (mod(nargin,2) == 0) | ...
        (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> fillinator:',...
        '    Invalid input arguments.',...
        '    Try "help fillinator".');
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
                'Error using ==> fillinator:',...
                '    Invalid property ', num2str(varargin{i}), ...
                ' is not (needs to be) a string.',...
                '    Try "help fillinator".');
            error(errorMessage);
        end
        fieldName = schnitzfield(varargin{i});
        p.(fieldName) = varargin{i+1};
    end
end

% Start out in the image directory to get names of existing images
eval(['cd ',p.imageDir]);
%Dy=dir([p.movieName,'-y-*.tif']);
%Dc=dir([p.movieName,'-c-*.tif']);
Dp=dir([p.movieName,'-p-*.tif']);

% get ACTUAL frames present (in case any are missing, starting at 0, etc)
if ~existfield(p,'fillinatorRange')
    numpos= findstr(Dp(1).name, '.tif')-3;
    [S,I] = sort({Dp.name}');
    frameNumStrings = char(S);
    p.fillinatorRange = str2num(frameNumStrings(:,numpos:numpos+2))';
end

fillstruc = [];
defaultlastframe = max(p.fillinatorRange);

notok = 1;
while notok,

    % Define the default cell number
    if ~existfield(p,'cell')
        p.cell = 1;
    end

    % Pompt user to enter the cell number, used in fillinator file name
    inputCellNum = input(['Cell number? [Default: ',num2str(p.cell),']']);
    if ~isempty(inputCellNum)
        p.cell = inputCellNum;
    end

    % Define the name of the default fillinator file
    if ~existfield(p,'fillinatorFile')
        p.fillinatorFile = [p.segmentationDir,'cell-',str3(p.cell),...
                            '-pixels.mat'];
    end
    
    % Prompt user to verify the name of the fillinator file to write
    disp(p.fillinatorFile);
    inputFilename = input(['Filename: [Default: (see above)] '],'s');
    if ~isempty(inputFilename),
        p.fillinatorFile = inputFilename;
    end;
    disp(['filename: ',p.fillinatorFile]);
    if exist(p.fillinatorFile),
        button = questdlg(['File exists, what do you want to do? ' ...
                           p.fillinatorFile],'Warning: file exists!',...
                           'Append','Overwrite','Cancel','Append');
        if upper(button(1))=='A',
            notok = 0;
            load(p.fillinatorFile,'fillstruc');
            defaultlastframe = -1+min([fillstruc.frame]);
        elseif upper(button(1))=='O',
            notok = 0;
        else
            notok = 1;
        end;
    else
        notok = 0;
    end;
    % verify we can write the fillinator file, in case filename screwed up
    if (notok == 0)
        fp = fopen(p.fillinatorFile,'w');
        if (fp == -1)
            disp(['Error: Could not write to file ' p.fillinatorFile]);
            notok = 1;
        end;
    end;
end;

last_frame = input(['Begin at? [Default: ',num2str(defaultlastframe),'] ']);
if isempty(last_frame),
    last_frame = defaultlastframe;
end;

% If user enters a last_frame that's not in fillinatorRange, let them know
if isempty(find(p.fillinatorRange==last_frame))
    error(['Frame ' num2str(last_frame) ' is not in fillinatorRange: ' ...
            num2str(p.fillinatorRange)]);
end
% Adjust default or argument fillinator range given user's input
p.fillinatorRange = sort(p.fillinatorRange(find(p.fillinatorRange<=last_frame)));

mainfig = figure;

set(mainfig,'WindowButtonMotionFcn',...
    'global mousepos;mousepos=get(gca,''CurrentPoint'');');

disp('After image has been processed...');
disp('    Left click to clear selection and repeat');
disp('    Press Space to accept, save, and continue');
disp('    Press ''S'' to skip and continue');
% disp('    Press ''a'' to increase boxsize and repeat');
% disp('    Press ''z'' to return boxsize to its initial (small) value and repeat');
% disp('    Press ''x'' to increase dilation and repeat');
% disp('    Press ''c'' to decrease dilation and repeat');
disp('    Press ''Y'' to use yfp as image');
disp('    Press ''C'' to use cfp as image');
disp('    Press ''P'' to use phase as image');
disp('    Press  ESC to reclick');
disp('    Press ''Q'' to quit');

maskcolors = 'YCP';
maxmask = 3;
minmask = 1;
minboxsize = 15;
maxboxsize = 200;
mask=3;
boxsize = 40;
defaultboxsize = 40;
frameoffset = 1;

numFrames = max(size(p.fillinatorRange));
% i is index into frame range, starting at end of range (highest frame number)
% loop will work through the fillinator range in reverse (high to low)
i = numFrames;
while i > 0
    frameNum = p.fillinatorRange(i);
    mynum = str3(frameNum);
    disp(['Processing frame ' mynum '...']);
   
    if mask>maxmask, mask=minmask; end;
    if mask<minmask, mask=maxmask; end;

    done=0;
    while ~done
        % supersize is because of the binning in the fluor images.  Ideally this
        % would be automatically determined instead of hard-coded
        if mask==1
            im = imread([p.imageDir p.movieName '-y-' mynum '.tif']);
            set(mainfig,'Name',['frame ' mynum ' YFP']);
            supersize = 2;
            invertit = 0;
        end
        if mask==2
            im = imread([p.imageDir p.movieName '-c-' mynum '.tif']);
            set(mainfig,'Name',['frame ' mynum ' CFP']);
            supersize = 2;
            invertit =0;
        end
        if mask==3
            im = imread([p.imageDir p.movieName '-p-' mynum '.tif']);
            supersize = 1;
            invertit = 1;
            set(mainfig,'Name',['frame ' mynum ' phase']);
        end

        disp('Click on a cell to flood it.');

        figure(mainfig);
        imshow(im,[]);
        title(['Click on a cell, or Right-click to go to menu']);
        [x,y,button] = ginput(1);
        x = round(x);
        y = round(y);
        subregion1 = max(y-boxsize/supersize,1):min(y+boxsize/supersize,size(im,1));
        subregion2 = max(x-boxsize/supersize,1):min(x+boxsize/supersize,size(im,2));
        subregion1 = round(subregion1);
        subregion2 = round(subregion2);
        newcoord1 = y - min(subregion1)+1;
        newcoord2 = x - min(subregion2)+1;

        box = im(subregion1,subregion2);
        if button == 1,
            [bw2] = progthreshsubtilis(box,0,supersize,invertit);
        else
            bw2 = zeros(size(box));
        end;

        bw2 = imopen(bw2,strel('disk',1));
        fullim = zeros(size(im));
        fullim(subregion1,subregion2) = bw2;
        imshow(makergb(im,fullim)); % this is the resulting image.
        title('<space>=OK+save, <esc>=redo, C,Y,P=switch image, S=skip, Q=quit, Bksp=previm, +/-/= adjust boxsize');
        % now wait for the command:

        didsomething = 0;
        while ~didsomething,
            w = waitforbuttonpress;
            while w==0,
                w = waitforbuttonpress;
            end;


            currChar = get(mainfig,'CurrentCharacter');
            disp('the char was')
            double(currChar)

            switch upper(currChar),
                case 'S'
                    done = 1;
                    didsomething =1;
                case 27 % ESC
                    didsomething = 1;
                case '+'
                    boxsize = round(boxsize*1.5);
                    if boxsize > maxboxsize,
                        boxsize = maxboxsize;
                    end;
                    boxsize,
                case '-'
                    boxsize = round(boxsize/1.5);
                    if boxsize < minboxsize,
                        boxsize = minboxsize;
                    end;
                    boxsize,
                case '=',
                    boxsize = defaultboxsize,
                case ' ',
                    % addandsave

                    if supersize > 1,
                        fullfullim = imresize(fullim,supersize);
                    else
                        fullfullim = fullim;
                    end;

                    thisone.source = maskcolors(mask)
                    [xx,yy] = find(fullfullim);
                    thisone.xpixels = xx;
                    thisone.ypixels = yy;
                    thisone.frame = frameNum;
                    if isempty(fillstruc),
                        fillstruc = thisone;
                    else
                        fillstruc(end+1) = thisone;
                    end;

                    f = find([fillstruc.frame]==thisone.frame);
                    if length(f)>1,
                        disp('there are duplicates, keeping most recent');
                        fillstruc(f(1:end-1))=[];
                    end;

                    if exist(p.fillinatorFile,'file')
                        copyfile(p.fillinatorFile,[p.fillinatorFile,'.bak']);
                    end
                    save(p.fillinatorFile, 'fillstruc');
                    didsomething=1;
                    done = 1;

                case 8, % backspace - go to back to LATER frame
                    if i < numFrames
                        i = i+2;     % will be decremented at end of loop
                    end

                    disp('prev image');
                    done=1;
                    didsomething=1;
                case 'Q'
                    button = questdlg ('Do you want to Quit?', 'Quit','No');
                    if upper(button(1))=='Y',
                        disp('Quitting');
                        return;
                    end;
                    didsomething=1;
                case {'C','Y','P'}
                    f = find(maskcolors==upper(currChar));
                    mask = f;
                    disp(['mask changed to: ',upper(currChar)]);
                    didsomething=1;
                case 'K'
                    keyboard;
                    didsomething=1;
            end;
        end;
    end;

    % decrement index into frame range and continue
    i = i-1;

end;  % while i > 0

