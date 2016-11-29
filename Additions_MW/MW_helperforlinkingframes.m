
function linklistschnitz =  MW_helperforlinkingframes(p, frameNumber)

%% Loading and resettings parameters
% Note that you standard 

% Which frame (n) are you interested in? (Will track n and n+1.)
%frameNumber=128;

% Load appropriate files for both n and n+1 (p is needed for this)
filename1 = [p.segmentationDir,p.movieName,'seg',str3(frameNumber),'.mat'];
load(filename1);
LcFrame1 = Lc;
filename2 = [p.segmentationDir,p.movieName,'seg',str3(frameNumber+1),'.mat'];
load(filename2);
LcFrame2 = Lc;

% Clear the linklist
linklist=[];

%% Linking frames
% Frame n in blue, frame n+1 in red.

% Create figure
h = figure(); clf; hold on;

% Start with labels on
showlabels=1;

while (any(LcFrame1(:)>0))

    disp('Note that matlab has variable editor in case you make mistake.');
    disp('Frame n is blue (first click), frame n+1 is red.');
    disp('Press ''l'' to link two cells or parent to daughter #2');
    disp('Press ''d'' to link parent to daughter #1');    
    disp('Press ''t'' to toggle labels off/on.');    
    disp('Press ''q'' to quit.');
    
    % Show them both
    MW_showtwoframeslabeled(LcFrame1, LcFrame2,h,showlabels)

    % Get input from user
    ct=waitforbuttonpress; % wait for buttonpress
    cc=get(h,'currentcharacter') % obtain value button % MW 2015/01
    
    % Terminate if 'q' pressed
    if cc=='q', break, end
    if cc=='t', showlabels=~showlabels; continue; end
    
    % Have user click such that 1st click is cell in frame 1 (n) and 2nd
    % click corresponds to cellno in frame 2 (n+1).
    click1 = round(ginput(1));
    cellno1 = LcFrame1(click1(2),click1(1));
    click2 = round(ginput(1));
    cellno2 = LcFrame2(click2(2),click2(1));

    % Update linklist
    linklist = [linklist; cellno1, cellno2]
    
    % Remove cell from image 
    if cc~='d' % Remove older cell if 'd' is not pressed
        LcFrame1(find(LcFrame1==cellno1)) = 0;
    end    
    LcFrame2(find(LcFrame2==cellno2)) = 0;
    
end

close(h)

%% Generate schnitz-compatible list

if ~isempty(linklist)
    linklistschnitz=MW_convertsimplelinklisttoschnitzlist(linklist,0)
    disp('Copy this list manually to corresponding file in data folder.')
else
    disp('Linklist was empty.');
end

disp('Bye.');






