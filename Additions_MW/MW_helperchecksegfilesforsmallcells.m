

%% Look for small cells in the frames
TRESHOLD = 20;

flag=0;
for frIdx = ourSettings.currentFrameRange
    
    %% 
    
    % Load segfile
    load([p.segmentationDir p.movieName 'seg' num2str(frIdx) '.mat'],'Lc');
    
    % Get all cell numbers for this frame
    allCellnos = unique(Lc(:))';
    
    % Loop over them to check size (brute force)
    for cellno = allCellnos
        if sum(sum(Lc==cellno)) < TRESHOLD
            disp(['cell ' num2str(cellno) ' in frame ' num2str(frIdx) ' is small']);
            flag=1;
        end
    end
    
    % Update user progress
    if mod(frIdx,100)==0
        disp(['Now at frame ' num2str(frIdx) ' of ' num2str(max(ourSettings.currentFrameRange))]);
    end
end
    
disp('Check done');    

%% How to fix it..

FRIDX=2152;
CELLNOsmall=155;

% load Lc
load([p.segmentationDir p.movieName 'seg' num2str(FRIDX) '.mat']);

% remove the culprit
Lc(Lc==CELLNOsmall)=0;

% renumber the Lc
allCellnos = unique(Lc(:));
allCellnosNonzero = allCellnos(allCellnos~=0)';
newNrs = 1:numel(allCellnosNonzero);
NewLc = zeros(size(Lc));
for idx = newNrs % goes from 1:nrCells 
    NewLc(Lc==allCellnosNonzero(idx)) = idx;
end
Lc=NewLc;
disp('Done re-indexing.');

% Loop cellnos to dubbel check size
allCellnos = unique(Lc(:))'; flag=0;
for cellno = allCellnos
    if sum(sum(Lc==cellno)) < TRESHOLD
        disp(['cell ' num2str(cellno) ' in frame ' num2str(frIdx) ' is small']);
        flag=1;
    end
end
if ~flag
    disp('No culprits found');
end

disp('Don''t forget to save!');

% and save
%save([p.segmentationDir p.movieName 'seg' num2str(FRIDX) '.mat'],'Lc','phsub','LNsub','rect','timestamp','phaseFullSize','tempsegcorrect');





