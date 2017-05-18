% CHANGED BY NOREEN
%
% DJK_getFluorShift takes original fluor images and segmentation files, and
% determines shift between fluor and phase image. Does this by trying
% different shifts and looking for maximum fluorescence within
% segmentation.
%
% OUTPUT
% 'optimalShift'  Best shift in case of getShift mode, otherwise empty.
%
% REQUIRED ARGUMENTS:
% 'p' 
%
% OPTIONAL ARGUMENTS:
% 'fluorcolor'    Fluor color that will be investigated. This can be
%                 'fluor1', 'fluor2', 'fluor3'. The colors associated with
%                 these expressions are stored in 'p'. If no color is
%                 given, 'fluor1'will be chosen.
% 'manualRange'   These frames will be treated
% 'DJK_saveDir'   Folder where tiff files are saved
%                 default: '\analysis\fluor\'
% 'maxShift'      Shifts the y image maximally this much for getShift 
%                 default: 7
% 'onScreen' = 0  automatically save and close images (default)
%              1  will ask to save and close images
% 'rescaleCorrection'  Rescales fluorescence image by an additional factor
%                 (apart from the factor due to different binning) and 
%                 recenters image (according to central image coordinates). In
%                 theory this should not be necessary, but for CFP images
%                 correction is needed (error in optics? etc) . Take around 
%                 XXXX for CFP.
%                 Default: =1

function optimalShift = DJK_getFluorShift_anycolor(p,varargin)

%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_getFluorShift_anycolor';

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
    % fluorcolor is special case, because the CONTENT of p.fluor1/2/3 must be written to
    % p.fluorcolor and not fluor1/2/3 literally. (NW 11/12/08)
    if strcmp(fieldName,'fluorcolor')==1 
        p.(fieldName)=p.(varargin{i+1});
        disp('--------------------------------------------------------------');
        disp(['Using ' varargin{i+1} ' (' p.(fieldName) ') as fluorescence color.']);
    else
        p.(fieldName) = varargin{i+1};
    end;
  end
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% If explicit manualRange is not given, take all segmentation files
if ~existfield(p,'manualRange')
  % Get directory of existing segmentation files 
  outprefix = [p.movieName 'seg'];
  D = dir([p.segmentationDir, outprefix, '*.mat']);
  [S,I] = sort({D.name}');
  D = D(I);
  numpos= findstr(D(1).name, '.mat')-3;
  
  segNameStrings = char(S);
  p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
end

% Choose fluorescence color if not specified
if ~isfield(p,'fluorcolor')
    p.fluorcolor=p.fluor1;
    disp('--------------------------------------------------------------');
    disp(['No fluorescence color specified. Will take fluor1 (' p.fluor1 ').']);
end

% If explicit DJK_saveDir is not given, define it
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'fluor_' p.fluorcolor filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end

% If explicit maxShift is not given, maxShift is 7
if ~existfield(p,'maxShift')
  p.maxShift = 7;
end

if ~existfield(p,'onScreen')
  p.onScreen = 0;
end

% If extra rescale Correction for Fluor Image does not exist, set to =1
% (probably all cases except for CFP)
if ~existfield(p,'rescaleCorrection')
    p.rescaleCorrection=1;
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% check if fluorcolor is non-existent ('none')
%--------------------------------------------------------------------------
if strcmp(p.fluorcolor,'none')==1
    disp('Fluorescence color is non-existent (set to ''none'').');
    optimalShift=[];
    return
else
    %--------------------------------------------------------------------------
    % Generate general variable names (e.g. yreg, ybinning, etc)
    %--------------------------------------------------------------------------
    reg=genvarname([p.fluorcolor 'reg']);
    binning=genvarname([p.fluorcolor 'binning']);
    back=genvarname([p.fluorcolor 'back']);
    gain=genvarname(['gain' p.fluorcolor]);
    expt=genvarname(['expt' p.fluorcolor]);
    
    
    %--------------------------------------------------------------------------
    % Open file to write results to
    %--------------------------------------------------------------------------
    frameShift = []; %holder for optimal shifts for each frame
    fid = fopen([p.DJK_saveDir  p.movieName '-getFluorShift.txt'],'wt');
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % Let know what is happening
    %--------------------------------------------------------------------------
    dispAndWrite(fid, ['-------------------------------------------------']);
    dispAndWrite(fid, ['Getting fluorShift from fluor images and seg files.']);
    dispAndWrite(fid, ['maxShift is ' num2str(p.maxShift)]);
    dispAndWrite(fid, ['extra rescale Correction is ' num2str(p.rescaleCorrection)]);
    dispAndWrite(fid, ['Only the (used) bestShifts directly from Fluorimages use the extra rescale correction!']);
    dispAndWrite(fid, ['-------------------------------------------------']);
    dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % LOOP OVER SEG FILES AND NORMALIZATION OF FLUOR
    %--------------------------------------------------------------------------
    % loop over frames
    for frameNum = p.manualRange

          % load complete fluor image (ALSO THIS WHEN p.TIFFonly)
          fluorname= [p.imageDir,p.movieName,'-',p.fluorcolor,'-',str3(frameNum),'.tif'];
          if exist(fluorname)==2
            fluorimage = imread(fluorname);
          else
            dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-' p.fluorcolor '-' str3(frameNum) '.tif in ' p.imageDir]);
            continue;
          end

          % load segmentation file
          filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];

          % The variables used to be cleared before reassignment but this is not possible with the
          % general names any more. Thus assignment of empty matrices. (NW 11/12/08)
          eval([reg '=[]; ' binning '=[]; ' back '=[]; ' gain '=[]; ' expt '=[]; ' ]);
          clear rect phaseFullSize Lc LNsub;
          load(filename);
          dispAndWrite(fid, [' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);
          if ~exist('Lc')      
            dispAndWrite(fid, ['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
            Lc = LNsub;
          end
          phaseCropSize = phaseFullSize; % seg file thinks phaseFullSize is full size, but might be crop size

          % still need to resize yimage
          eval(['fluorimage = imresize_old(fluorimage,' binning ',''nearest'');']);
          
         
          %------------------------------------------------------------------------
          % perform manual rescale Correction (NW_2012-02-27)
          %------------------------------------------------------------------------
          if (length(p.rescaleCorrection==1) & p.rescaleCorrection~=1) | length(p.rescaleCorrection)>1
              centerfluor=NW_rescalecenterimage_affine(fluorimage,p.rescaleCorrection);
              fluorimage=centerfluor;
              clear centerfluor;              
          end
          %------------------------------------------------------------------------          
 
                    
          %------------------------------------------------------------------------
          % Check Lc & [y/c/...]reg from segmentation file
          %------------------------------------------------------------------------
          % get subimages
          eval(['dblfluor = double(' reg ');']);
          seg = +(Lc(1+ p.maxShift:end- p.maxShift, 1+ p.maxShift:end- p.maxShift) > 0);
          % try all possible translations
          for di = - p.maxShift: p.maxShift
            for dj = - p.maxShift: p.maxShift
              fluor_shifted = dblfluor( p.maxShift+di+1:end- p.maxShift+di,  p.maxShift+dj+1:end- p.maxShift+dj);
              tot_fluor(di+ p.maxShift+1,dj+ p.maxShift+1) = sum(fluor_shifted(seg > 0));
            end
          end
          
          % best translation is the one with largest fluorescence within seg
          [shift_y, shift_x] = find(tot_fluor == max2(tot_fluor)); % find returns [row, col]
          best_shift(1) = - (mean(shift_x) -  p.maxShift - 1);
          best_shift(2) = - (mean(shift_y) -  p.maxShift - 1);
          %------------------------------------------------------------------------

          %------------------------------------------------------------------------
          % Check Lc & loaded fluorimage (these data for best shift will be
          % output of fct)
          %------------------------------------------------------------------------
          % get subimages
          dblfluor=double(fluorimage);
          L = zeros(phaseCropSize);
          L((rect(1):rect(3)), (rect(2):rect(4))) = Lc;
          seg = +(L(1+ p.maxShift:end- p.maxShift, 1+ p.maxShift:end- p.maxShift) > 0);
          % try all possible translations
          for di = - p.maxShift: p.maxShift
            for dj = - p.maxShift: p.maxShift
              fluor_shifted = dblfluor( p.maxShift+di+1:end- p.maxShift+di,  p.maxShift+dj+1:end- p.maxShift+dj);
              tot_fluor(di+ p.maxShift+1,dj+ p.maxShift+1) = sum(fluor_shifted(seg > 0));
            end
          end
          % best translation is the one with largest fluorescence within seg
          [shift_y, shift_x] = find(tot_fluor == max2(tot_fluor));
          best_shift(3) = - (mean(shift_x) -  p.maxShift - 1);
          best_shift(4) = - (mean(shift_y) -  p.maxShift - 1);

          % record and display translation
          frameShift(size(frameShift,1)+1,:) = [frameNum best_shift(1) best_shift(2) best_shift(3) best_shift(4)];
          dispAndWrite(fid, ['       -> best shift of fluor is [' num2str(best_shift(1)) ',' num2str(best_shift(2)) '] and [' num2str(best_shift(3)) ',' num2str(best_shift(4)) ']']);
          %------------------------------------------------------------------------

    end 
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % AVERAGE OUTPUT OF SHIFT CHECK
    %--------------------------------------------------------------------------
    if ~(size(frameShift)==[0,0])
        
        average = round(mean(frameShift));
        dispAndWrite(fid, ['-------------------------------------------------']);
        if numel(average)>=5
            dispAndWrite(fid, [' * Best average shifts are  [' num2str(average(2)) ',' num2str(average(3)) '] and  [' num2str(average(4)) ',' num2str(average(5)) ']']);
        else
            warning('average parameter has too little information');
        end

        % Show histograms
        scrsz = get(0, 'ScreenSize');
        fig1 = figure('Position', [151 scrsz(4)-200 scrsz(3)-130 scrsz(4)-200], 'visible','off');
        hist(frameShift(:,4));
        xlabel('Optimal shifts for x');
        xlim([-p.maxShift-1 p.maxShift+1]);
        hold off;

        fig2 = figure('Position', [151 scrsz(4)-200 scrsz(3)-130 scrsz(4)-200], 'visible','off');
        hist(frameShift(:,5));
        xlabel('Optimal shifts for y');
        xlim([-p.maxShift-1 p.maxShift+1]);

        % Ask to save the figure
        if p.onScreen
          set(fig1,'visible','on');
          set(fig2,'visible','on');
          saveFigInput = questdlg('Save Figures?','Save Figures?','Yes','Yes and Close','No','Yes');
          pause(0.2);
        else
          saveFigInput='Yes and Close';
        end

        if (upper(saveFigInput(1))=='Y')
        %   saveas(fig1,[p.DJK_saveDir p.movieName '-histogram_optimal_shift_x.fig']);
        %   saveas(fig2,[p.DJK_saveDir p.movieName '-histogram_optimal_shift_y.fig']);
          saveSameSize(fig1,'file',[p.DJK_saveDir p.movieName '-histogram_optimal_shift_x.png'], 'format', 'png');
          saveSameSize(fig2,'file',[p.DJK_saveDir p.movieName '-histogram_optimal_shift_y.png'], 'format', 'png');
          if (strcmp(saveFigInput,'Yes and Close'))
            close(fig1);close(fig2);
            pause(0.2);
          end
          dispAndWrite(fid, [' * Saved histogram of optimal shift for x in ' p.movieName '-histogram_optimal_shift_x.png']);
          dispAndWrite(fid, [' * Saved histogram of optimal shift for y in ' p.movieName '-histogram_optimal_shift_y.png']);
        end
    else
        warning('Optimal shift not found, size(frameShift)==[0,0]. Setting shift to [0,0]. TODO.');
        disp('Waiting five seconds before continuing..');
        % I'm not sure what goes wrong when this happens. But I should look
        % into how this function works. For now I choose the boilerplate
        % solition, which is ignore if the function spits out nonsense, and
        % just perform no shift..
        average=zeros(1,5);
        pause(5);
    end
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % Close file to write results to
    %--------------------------------------------------------------------------
    dispAndWrite(fid, ['-------------------------------------------------']);
    fclose(fid);
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % OUTPUT
    %--------------------------------------------------------------------------
    optimalShift = [average(4) average(5)];
    %--------------------------------------------------------------------------


end
