function DJK_analyzeFluorBackground_anycolor(p,varargin)


% Function has been marked as obsolete.

warning('This function has been marked as obsolete');




% function DJK_analyzeFluorBackground_anycolor(p,varargin)
% CHANGED BY NOREEN
%
% loads the extra rescaled corrected fluorescence data from .mat file resp.
% the not corrected data from original fluorimage (NW 2012-02-27) (image directory) 
% resp corrected fluorimages from the directories c2, c3, c4 etc
%
% DJK_analyzeFluorBackground_anycolor
%
% OUTPUT
%
% REQUIRED ARGUMENTS:
%
% OPTIONAL ARGUMENTS:
% 'fluorcolor'    Fluor color that will be investigated. This can be
%                 'fluor1', 'fluor2', 'fluor3'. The colors associated with
%                 these expressions are stored in 'p'. If no color is
%                 given, 'fluor1'will be chosen.
% 'manualRange'   These frames will be treated
%
% 'DJK_saveDir'   Folder where tiff files are saved
%                 (default: '\analysis\fluor_[color]\')
%
% 'TIFFonly' = 1  does nothing with seg files, only creates a [color]2 (e.g. y2) TIFF 
%                 manualRange needs to be given, else might give
%                 error when loading seg file
%
% 'onScreen' = 1  will not automatically save and close images, but ask
% 'rescaleCorrection'  Rescales fluorescence image by an additional factor
%                 (apart from the factor due to different binning) and 
%                 recenters image (according to central image coordinates). In
%                 theory this should not be necessary, but for CFP images
%                 correction is needed (error in optics? etc) . Take around 
%                 XXXX for CFP.
%                 Default: =1
% 'minimalMode'       =1: only calculates values, that are actually used
%                     later (Y5xxx,Y6_mean etc).
%                     =0: calculates everything (default: =0)  (NW 2012/04)
%


%--------------------------------------------------------------------------
% Input error checking
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_analyzeFluorBackground_anycolor';

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
% make sure every directory field has a trailing filesep
if (p.DJK_saveDir(length(p.DJK_saveDir)) ~= filesep)
  p.DJK_saveDir = [p.DJK_saveDir filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end

% If explicit TIFFonly is not given, TIFFonly is false (0)
if ~existfield(p,'TIFFonly')
  p.TIFFonly = 0;
end

% Just in case it has not been set yet
if ~existfield(p,'onScreen')
  p.onScreen = 0;
end

% If extra rescale Correction for Fluor Image does not exist, set to =1
% (probably all cases except for CFP)
if ~existfield(p,'rescaleCorrection')
    p.rescaleCorrection=1;
end
if ~existfield(p,'minimalMode') % (NW 2012/04)
  p.minimalMode = 0;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Open file to write results to
%--------------------------------------------------------------------------
fid = fopen([p.DJK_saveDir p.movieName '-fluorBackground.txt'],'wt');
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Let know what is happening
%--------------------------------------------------------------------------
dispAndWrite(fid, ['-------------------------------------------------']);
dispAndWrite(fid, ['Analyzing ' num2str(length(p.manualRange)) ' frames from ' num2str(p.manualRange(1)) ' to ' num2str(p.manualRange(end))]);
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
    % Generate general variable names (e.g. all_rback, all_rbackNew3, etc)
    %--------------------------------------------------------------------------
    all_back=genvarname(['all_' p.fluorcolor 'back']);           % e.g. all_rback
    all_back2=genvarname(['all_' p.fluorcolor 'back2']);         % e.g. all_rback2
    all_back2Alt=genvarname(['all_' p.fluorcolor 'back2Alt']);   % e.g. all_rback2Alt
    all_back3=genvarname(['all_' p.fluorcolor 'back3']);         % e.g. all_rback3
    all_backAlt=genvarname(['all_' p.fluorcolor 'backAlt']);     % e.g. all_rbackAlt
    all_backNew=genvarname(['all_' p.fluorcolor 'backNew']);     % e.g. all_rbackNew
    all_backNew2=genvarname(['all_' p.fluorcolor 'backNew2']);   % e.g. all_rbackNew2
    all_backNew3=genvarname(['all_' p.fluorcolor 'backNew3']);   % e.g. all_rbackNew3
    reg=genvarname([p.fluorcolor 'reg']);
    binning=genvarname([p.fluorcolor 'binning']);
    back=genvarname([p.fluorcolor 'back']);
    gain=genvarname(['gain' p.fluorcolor]);
    expt=genvarname(['expt' p.fluorcolor]);
    back2=genvarname([p.fluorcolor 'back2']);
    back3=genvarname([p.fluorcolor 'back3']);
    backAlt=genvarname([p.fluorcolor 'backAlt']);
    back2Alt=genvarname([p.fluorcolor 'back2Alt']);
    backNew=genvarname([p.fluorcolor 'backNew']);
    backNew2=genvarname([p.fluorcolor 'backNew2']);
    backNew3=genvarname([p.fluorcolor 'backNew3']);
    
    %--------------------------------------------------------------------------
    % LOOP OVER SEG FILES AND GET BACKGROUND DATA
    %--------------------------------------------------------------------------
    all_frames = [];
    eval([all_back '= [];']);
    eval([all_back2 '= [];']);
    eval([all_back2Alt '= [];']);
    eval([all_back3 '= [];']);
    eval([all_backAlt '= [];']);
    eval([all_backNew '= [];']);
    eval([all_backNew2 '= [];']);
    eval([all_backNew3 '= [];']);


    % loop over frames
    count = 0;
    for frameNum = p.manualRange
      %------------------------------------------------------------------------
      % LOAD DATA
      %------------------------------------------------------------------------
      % load complete fluor image (not extra rescale corrected)
      fluorname= [p.imageDir,p.movieName,'-' p.fluorcolor '-',str3(frameNum),'.tif'];
      if exist(fluorname)==2
        fluorimage = imread(fluorname);
      else
        disp([' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-' p.fluorcolor '-' str3(frameNum) '.tif in ' p.imageDir]);
        continue;
      end

      % load complete fluor image (already extra rescale corrected)
      fluor2name= [p.analysisDir, 'fluor_', p.fluorcolor filesep p.fluorcolor '2' filesep, p.movieName ,'-' p.fluorcolor '2-',str3(frameNum),'.tif'];
      if exist(fluor2name)==2
        fluor2image = imread(fluor2name);
      else
        dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-' p.fluorcolor '2-' str3(frameNum) '.tif in ' p.analysisDir, 'fluor_', p.fluorcolor filesep p.fluorcolor '2' filesep]);
        continue;
      end

      if ~p.TIFFonly %(already extra rescale corrected)
            fluor3name= [p.analysisDir, 'fluor_', p.fluorcolor filesep p.fluorcolor '3' filesep, p.movieName ,'-' p.fluorcolor '3-',str3(frameNum),'.tif'];
            if exist(fluor3name)==2
               fluor3image = imread(fluor3name);
            else
               dispAndWrite(fid, [' * ' str3(frameNum) ' -> Cannot open : ' p.movieName '-' p.fluorcolor '3-' str3(frameNum) '.tif in ' p.analysisDir, 'fluor_', p.fluorcolor filesep p.fluorcolor '3' filesep]);
               continue;
            end
          
            % load segmentation file
            filename = [p.segmentationDir, p.movieName, 'seg', str3(frameNum)];
            % The variables used to be cleared before reassignment but this is not possible with the
            % general names any more. Thus assignment of empty matrices. (NW 11/12/08)
            % NOTE: It seems as if the variables could indeed be cleared. (NW 11/12/13)
            eval([reg '=[]; ' binning '=[]; ' back '=[]; ' gain '=[]; ' expt '=[]; ' ]);
            clear rect phaseFullSize Lc LNsub;
            load(filename);
            dispAndWrite(fid, [' * ' str3(frameNum) ' -> loaded ' p.movieName 'seg' str3(frameNum) ' in ' p.segmentationDir]);

            % load fluor file
            filename = [p.tracksDir p.movieName 'Fluor_' p.fluorcolor '_' str3(frameNum)];
            load(filename);
            dispAndWrite(fid, ['       -> loaded ' p.movieName 'Fluor_' p.fluorcolor '_' str3(frameNum) ' in ' p.tracksDir]);

            % get segmented part correct
            if ~exist('Lc')      
              dispAndWrite(fid, ['       ->  segmentations has not been corrected -> will use LNsub in stead of Lc !!!']);
              Lc = LNsub;
            end
            segcrop = +(Lc > 0);
            seg = zeros(phaseCropSize);
            seg(rect(1):rect(3),rect(2):rect(4)) = segcrop;
      else
            % Get sizes correct
            phaseCropSize = phaseFullSize;
            seg = zeros(phaseCropSize);
      end

      % still need to resize [color]image
      eval(['fluorimage = imresize_old(fluorimage,' binning ',''nearest'');']);
      %------------------------------------------------------------------------
      % perform manual rescale Correction of fluorimage(NW 2012-02-27)
      %------------------------------------------------------------------------
      if (length(p.rescaleCorrection==1) & p.rescaleCorrection~=1) | length(p.rescaleCorrection)>1
          centerfluor=NW_rescalecenterimage_affine(fluorimage,p.rescaleCorrection);
          fluorimage=centerfluor;
          clear centerfluor;              
      end
      %------------------------------------------------------------------------    
      % fluor2image was corrected in DJK_correctFluorImage_anycolor for
      % extra rescaling
      eval(['fluor2image = imresize_old(fluor2image,' binning ',''nearest'');']);
      %------------------------------------------------------------------------
      
      
      %------------------------------------------------------------------------
      % COLLECT BACKGROUND DATA
      %------------------------------------------------------------------------
      all_frames = [all_frames frameNum];
      eval([all_back '= [' all_back ' ' back '];']);
      if p.minimalMode~=1
          eval([all_back2 '= [' all_back2 ' ' back2 '];']);
          eval([all_back2Alt  '= [' all_back2Alt ' ' back2Alt '];']);
          eval([ all_backAlt '= [' all_backAlt ' ' backAlt '];']);
      end
      eval([all_back3 '= [' all_back3 ' ' back3 '];']);
      
      %------------------------------------------------------------------------


      if p.minimalMode~=1
          %------------------------------------------------------------------------
          % GET NEW BACKGROUND DATA 1
          %------------------------------------------------------------------------
          seg_dilate = imdilate(seg,strel('disk',30));
          [rows, cols] = find(seg_dilate>0);
          row_min = min(rows);
          row_max = max(rows);
          col_min = min(cols);
          col_max = max(cols);
          rect_dilate = [row_min col_min row_max col_max];

          background1 = zeros(phaseCropSize);
          background1(row_min:row_max,col_min:col_max) = 1;
          background1(find(seg_dilate>0)) = 0;
          idx = find(background1>0);

          eval([ backNew ' = median(  fluor3image(idx));']);
          eval([all_backNew '= [' all_backNew ' ' backNew '];']);
          %------------------------------------------------------------------------


          %------------------------------------------------------------------------
          % GET NEW BACKGROUND DATA 2 & 3
          %------------------------------------------------------------------------
          % first get boundary of microcolonie
          seg_filled = imfill(imclose(seg, strel('disk',1)));
          seg_filled_perim = bwperim(seg_filled);
          [perim_y, perim_x] = find(seg_filled_perim>0);

          % get min distance for each pixel away from seg
          [nonseg_y, nonseg_x] = find(seg_filled==0);
          nonseg_idx = find(seg_filled==0);

          %disp('starting calc optimized'); t=cputime; 
          for i = 1:length(perim_y)
            dist_squared = (nonseg_y-perim_y(i)).^2 + (nonseg_x-perim_x(i)).^2;
            if i==1, min_dist_squared = dist_squared;
            else, min_dist_squared = min(min_dist_squared, dist_squared);
               % if i==100, disp(min_dist_squared(1:100)); end % blubb
            end
          end
          %disp(['finished calc optimized in ' num2str(cputime-t)]); pause(0.2);
          % keyboard;

          pixels_squared_per_micro = (1/p.micronsPerPixel)^2;

          % 0 tot 1 um etc
          for i = 1:7
            pixels2 = fluor3image( nonseg_idx(find(min_dist_squared>=(i-1)*pixels_squared_per_micro & min_dist_squared<=i*pixels_squared_per_micro)) );
            pixels3 = fluor2image( nonseg_idx(find(min_dist_squared>=(i-1)*pixels_squared_per_micro & min_dist_squared<=i*pixels_squared_per_micro)) );
            if length(pixels2) > 0
                eval([backNew2 '(i) = median(pixels2);']);
                eval([backNew3 '(i) = median(pixels3);']);
            else
                eval([backNew2 '(i) = 0;'])
                eval([backNew3 '(i) = 0;'])
            end
          end
          count = count + 1;
          eval([all_backNew2 '(:,count) = [' backNew2 '(1); ' backNew2 '(2); ' backNew2 '(3); ' backNew2 '(4); ' backNew2 '(5); ' backNew2 '(6); ' backNew2 '(7); ];']);
          eval([all_backNew3 '(:,count) = [' backNew3 '(1); ' backNew3 '(2); ' backNew3 '(3); ' backNew3 '(4); ' backNew3 '(5); ' backNew3 '(6); ' backNew3 '(7); ];']);
      end %end p.minimalMode~=1
          %------------------------------------------------------------------------
    end %end loop frames
        %--------------------------------------------------------------------------

        %--------------------------------------------------------------------------
        % PRINT BACKGROUND DATA
        %--------------------------------------------------------------------------
        dispAndWrite(fid, ['-------------------------------------------------']);
        eval(['dispAndWrite(fid, [''frame ' back ' ' back2 ' ' back2Alt ' ' back3 ' ' backAlt ' ' backNew ' 0-1 1-2 2-3 3-4 4-5 5-6 6-7 '']);']);
        for i = 1:length(all_frames)
            if p.minimalMode~=1
                eval(['dispAndWrite(fid, [ str3(all_frames(i)) '' ''  str3(' all_back '(i)) '' '' str3(' all_back2 '(i)) '' ''  str3(' all_back2Alt '(i)) '' '' str3(' all_back3 '(i)) '' ''  str3(' all_backAlt '(i)) '' ''  str3(' all_backNew '(i)) '' ''  str3(' all_backNew2 '(1,i)) '' ''  str3(' all_backNew2 '(2,i)) '' '' str3(' all_backNew2 '(3,i)) '' '' str3(' all_backNew2 '(4,i)) '' '' str3(' all_backNew2 '(5,i)) '' ''  str3(' all_backNew2 '(6,i)) '' '' str3(' all_backNew2 '(7,i)) '' '' str3(' all_backNew3 '(1,i)) '' '' str3(' all_backNew3 '(2,i)) '' '' str3(' all_backNew3 '(3,i)) '' '' str3(' all_backNew3 '(4,i)) '' '' str3(' all_backNew3 '(5,i)) '' ''  str3(' all_backNew3 '(6,i)) '' '' str3(' all_backNew3 '(7,i)) '' '' ]);']);
            end
        end
        %--------------------------------------------------------------------------


        %--------------------------------------------------------------------------
        % NORMALIZE FLUOR DATA SAME WAY AS FOR SCHNITZES
        %--------------------------------------------------------------------------
        if p.minimalMode~=1
            eval([all_back2 '= normalize( ' all_back2 ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
            eval([all_back2Alt '= normalize( ' all_back2Alt ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
            eval([all_backAlt '= normalize( ' all_backAlt ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
            eval([all_backNew '= normalize( ' all_backNew ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
            eval([all_backNew2 '= normalize( ' all_backNew2 ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
            eval([all_backNew3 '= normalize( ' all_backNew3 ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
        end 
    eval([all_back '= normalize( ' all_back ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
    eval([all_back3 '= normalize( ' all_back3 ', 0, ' expt ', ' binning ', p.micronsPerPixel, 0);']);
       
    
    %--------------------------------------------------------------------------

    if p.minimalMode~=1
        %--------------------------------------------------------------------------
        % Show Plot
        %--------------------------------------------------------------------------
        % make figure
        scrsz = get(0, 'ScreenSize');
        fig111 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
        hold on;

        subplot(2,1,1);
        eval(['plot(all_frames, ' all_back ', ''b-'', all_frames, ' all_backAlt ',''k:''); ']);
        legend(back,backAlt,-1);
        %ylim([10 15]); blubb

        subplot(2,1,2);
        eval(['plot(all_frames, ' all_back2 ', ''b-'', all_frames, ' all_back2Alt ',''k:'', all_frames, ' all_back3 ',''b--'',all_frames, ' all_backNew ',''r-.'');']); 
        legend(back2,back2Alt,back3,backNew,-1);
        %ylim([0 5]); blubb
        xlabel('Frames');

        fig112 = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
        hold on;

        subplot(2,1,1);
        eval(['plot( all_frames, ' all_backNew2 '(1,:), ''g-'', all_frames, ' all_backNew2 '(2,:), ''r-'', all_frames, ' all_backNew2 '(3,:), ''b-'', all_frames, ' all_backNew2 '(4,:), ''g--'', all_frames, ' all_backNew2 '(5,:), ''r--'', all_frames, ' all_backNew2 '(6,:), ''b--'', all_frames, ' all_backNew2 '(7,:), ''g-.''); ']);
        legend([backNew2 ' 0 to 1 um'],[backNew2 ' 1 to 2 um'],[backNew2 ' 2 to 3 um'],[backNew2 ' 3 to 4 um'],[backNew2 ' 4 to 5 um'],[backNew2 ' 5 to 6 um'],[backNew2 ' 6 to 7 um'],-1);
        % ylim([0 5]); blubb

        subplot(2,1,2);
       eval(['plot( all_frames, ' all_backNew3 '(1,:), ''g-'', all_frames, ' all_backNew3 '(2,:), ''r-'',  all_frames, ' all_backNew3 '(3,:), ''b-'', all_frames, ' all_backNew3 '(4,:), ''g--'', all_frames, ' all_backNew3 '(5,:), ''r--'', all_frames, ' all_backNew3 '(6,:), ''b--'',  all_frames, ' all_backNew3 '(7,:), ''g-.''); ']);
        legend([backNew3 ' 0 to 1 um'],[backNew3 ' 1 to 2 um'],[backNew3 ' 2 to 3 um'],[backNew3 ' 3 to 4 um'],[backNew3 ' 4 to 5 um'],[backNew3 ' 5 to 6 um'],[backNew3 ' 6 to 7 um'],-1);
        % ylim([0 5]); blubb
        xlabel('Frames');

        % Ask to save the figure
        if p.onScreen
          saveFigInput = questdlg('Save Figures?','Save Figures?','Yes','Yes and Close','No','Yes');
          pause(0.2);
        else
          saveFigInput='Yes and Close';
        end

        if (upper(saveFigInput(1))=='Y')
          saveas(fig111,[p.DJK_saveDir p.movieName '-fluor_' p.fluorcolor '_Background1.fig']);
          saveSameSize(fig111,'file',[p.DJK_saveDir p.movieName '-fluor_' p.fluorcolor '_Background1.png'], 'format', 'png');
          saveas(fig112,[p.DJK_saveDir p.movieName '-fluor_' p.fluorcolor '_Background2.fig']);
          saveSameSize(fig112,'file',[p.DJK_saveDir p.movieName '-fluor_' p.fluorcolor '_Background2.png'], 'format', 'png');
          if (strcmp(saveFigInput,'Yes and Close'))
            close(fig111); close(fig112);
            pause(0.2);
          end
          dispAndWrite(fid, ['-------------------------------------------------']);
          dispAndWrite(fid, [' * Saved plot background fluor (' p.fluorcolor ') over time in ' p.movieName '-fluor_' p.fluorcolor '_Background1.png']);
          dispAndWrite(fid, [' * Saved plot background fluor (' p.fluorcolor ') over time in ' p.movieName '-fluor_' p.fluorcolor '_Background2.png']);
        end
        %--------------------------------------------------------------------------
    end % end p.minimalMode~=1


    %--------------------------------------------------------------------------
    % Close file to write results to
    %--------------------------------------------------------------------------
    dispAndWrite(fid, ['-------------------------------------------------']);
    if fid~=0
      fclose(fid);
    end
    %--------------------------------------------------------------------------
    
    warning('This function has been marked as obsolete');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize performs new normalization
function fluorregNorm = normalize( fluorreg, ...            % Fluor image in camera units
                               fluorback, ...           % background in camera units
                               exptfluor, ...           % exposure time in ms
                               fluorbinning, ...        % camera bin size (2)
                               micronsPerPixel, ... % micron per pixel (0.04065)
                               autoFluor)           % autofluorescence of cells
if (exptfluor=='emptr')
  disp('Warning! no exposure time specified for fluor image. Arbitrarily set to 1000');
  exptfluor = 1000; 
end
fluorregNorm=double(fluorreg)-double(fluorback);
fluorregNorm=fluorregNorm/(exptfluor * fluorbinning^2 * micronsPerPixel^2);
fluorregNorm=fluorregNorm-autoFluor; % in units of molecules?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


