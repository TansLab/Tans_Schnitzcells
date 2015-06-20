% TODO MW
% ===
%
% This function is depricated and should not be used!
% 
% It does contain some code to accomodate a 7-degree polynomial fit 
% for the length, but rather not use different funtions for different
% utilities. 
% TODO: recode DJK_addToSchnitzes_length such that it will contain an
% option to use 7-degree polynomial.
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% This file is a quick and dirty merge of Rutgers code into the original
% DJK_addToSchnitzes_length.m. This file should be thoroughly checked
% whether it contains no bugs.

% DJK_addToSchnitzes_length loads schnitzcells, calculates some better
% lengths of cells, and saves them to schnitzcells again.
%
% DJK_compileSchnitzImproved must have been run before 
%
% Option to plot is available. Will plot data for every 100tht schnitz
% only, for more use schnitzNum array.
%
% OUTPUT
% 'p'   
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'p'
%
% OPTIONAL ARGUMENTS:
% 'schnitzName'       
% 'micronsPerPixel'   default = 0.04065 (1.5x magnification of Killer Mike)
% 'onScreen' = 0      automatically save and close images (default only
%                     schnitz 1, see p.schnitzNum parameter below)
%              1      will ask to save and close images
%              2      will not show or save images (default)
% 'schnitzNum'        Array with schnitz numbers whose data will be plotted
%                     default: 1st schnitz only
% 'DJK_saveDir'       directory where images will be saved. 
%                     default: [p.analysisDir 'schnitzLength\']
%
% p.schnitzNum        If onScreen = 0 or 1, one can use p.schnitzNum=[a..b] 
%                     to give range of which Schnitzes to save as figure. 
%                     Keyword p.schnitzNum='all' will output all of them 
%                     (if onscreen = 0 or 1).
%
%

function [p,schnitzcells] = MW_addToSchnitzes_length(p,varargin) 

disp('*ERROR*: This function should not be used; use DJK_addToSchnitzes_length instead');
error('*ERROR*: This function should not be used; use DJK_addToSchnitzes_length instead');
return % comment to forcefully use

%--------------------------------------------------------------------------
% Input error checking and parsing
%--------------------------------------------------------------------------
% Settings
numRequiredArgs = 1; functionName = 'DJK_addToSchnitzes_length';

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
% Set default parameter values if they don't exist yet
if ~existfield(p,'schnitzName')
  p.schnitzName = [p.tracksDir,p.movieName,'-Schnitz.mat'];
end
if ~existfield(p,'micronsPerPixel')
  p.micronsPerPixel = 0.04065;
end
% If onScreen, nothing is saved to disc automatically
if ~existfield(p,'onScreen')
  p.onScreen = 2;
end
% Save directory
if ~existfield(p,'DJK_saveDir')
  p.DJK_saveDir = [p.analysisDir 'schnitzLength' filesep];
end
% Make sure that DJK_saveDir directory exists
if exist(p.DJK_saveDir)~=7
  [status,msg,id] = mymkdir([p.DJK_saveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' p.DJK_saveDir ' : ' msg]);
    return;
  end
end
%--------------------------------------------------------------------------
  

%--------------------------------------------------------------------------
% Load Schnitzcells
%--------------------------------------------------------------------------
% Loading an existing schnitz file which contains the image-derived fields
if exist(p.schnitzName)~=2
  error(['Could not read ' p.schnitzName ' , which is required for quick mode']);
end
load(p.schnitzName);
disp(['Load from ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Overwrite any schnitzcells parameters/defaults given optional fields/values
%--------------------------------------------------------------------------
% Reintroducing option to plot all fits - MW 23/3/2014
if existfield(p,'schnitzNum')  
    if p.schnitzNum == 'all'
        p.schnitzNum = [1:length(schnitzcells)];
        % p.schnitzNum = [1:100:length(schnitzcells)]; %; %blubb
    end   
else
 p.schnitzNum = [1];
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% CREATION OF lincellnum
%--------------------------------------------------------------------------
% Will use a lincellnum structure, where for only used frames (not necessarly
% starting with image 001 (frame 2), the link between schnitznum and cellno
% is saved: lincellnum {a} (b) returns schnitznum of cellno b in frame a

% Get trackRange (frames that will be extracted). 
trackRange = sort(unique([schnitzcells.frame_nrs])); % MW 2014/06/11 removal N+1 bug

% initialize lincellnum to have zero'd arrays for each frame 
lincellnum = {};
for lincellnumIndex = 1:length(trackRange)
  segdata = load([p.segmentationDir,p.movieName,'seg',str3(trackRange(lincellnumIndex))]);
  if ~isfield(segdata,'Lc');
    segdata.Lc = segdata.LNsub;
  end
  lincellnum{lincellnumIndex} = zeros([1 max2(segdata.Lc)]);
end

% Now, step through each schnitz, store schnitz number in lincellnum
for schnitznum = 1:length(schnitzcells)
  s = schnitzcells(schnitznum);
  for age = 1:length(s.frame_nrs)
    framenum = s.frame_nrs(age); % MW 2014/06/11 removal N+1 bug

    % look up that frame number's index within trackRange
    lincellnumIndex = find(trackRange==framenum);
    cellnum  = s.cellno(age);
    lincellnum{lincellnumIndex}(cellnum) = schnitznum;
  end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER FRAMES IN TRACKRANGE
%--------------------------------------------------------------------------
for i = 1:length(trackRange)
  currFrameNum = trackRange(i);

  % i is the frame index into lincellnum and trackRange. 
  % trackRange(i) gives the actual frame number that corresponds 
  % to a given lincellnum{i}.
  % Also, recall that a schnitz' frames array is 1-based per 
  % makeschnitz hack to conform to schnitzedit expectations

  % empty previous data
  clear Lc phsub;
  timestamp = 0;

  % load segmented image for this frameNum
  name = [p.segmentationDir,p.movieName,'seg',str3(currFrameNum)];
  disp([' Processing frame ', str3(currFrameNum), ' as nr ', str3(i), ' in range']);
  load([name]); % including variables: LNsub, Lc, phsub, timestamp, rect, yback, yreg, yshift, expty, gainy    

  %------------------------------------------------------------------------
  % LOOP OVER EACH SCHNITZ THAT EXISTS DURING THIS FRAME, AND UPDATE IT
  %------------------------------------------------------------------------
  schnitzesForFrame = lincellnum{i}; % [schnitznum for cellno1, schnitznum for cellno2, etc]
  nonZeroSchnitzes = schnitzesForFrame(schnitzesForFrame~=0); % with correction you sometimes end up with unexisting schnitzes (0)
  
  
  for s = nonZeroSchnitzes
    % figure out index within this schnitz' age-based arrays
    age = find((schnitzcells(s).frame_nrs) == currFrameNum); % MW 2014/06/11 removal N+1 bug
    if isempty(age)
      error(['lincellnum says schnitz num ' s 'exists in frame ' currFrameNum ', but that frame can''t be found in the schnitz' frames array']);
    end

    % get cell number in segmented image for current schnitz & frame
    cellnum = schnitzcells(s).cellno(age);

    % get previous stuff
    phi = schnitzcells(s).phi(age);
    fitCoef2 = schnitzcells(s).fitCoef2(age,:);
    fitCoef3 = schnitzcells(s).fitCoef3(age,:);
    fitCoef4 = schnitzcells(s).fitCoef4(age,:); % RUTGER ADDED
   
    % HIGHER ORDER FIT
    fitCoef5 = schnitzcells(s).fitCoef5(age,:); % RUTGER ADDED
    fitCoef6 = schnitzcells(s).fitCoef6(age,:); %
    fitCoef7 = schnitzcells(s).fitCoef7(age,:); %
 %  func_3rd = @(x) x.^3 .* fitCoef3(1) + x.^2 .* fitCoef3(2) + x .* fitCoef3(3) + fitCoef3(4);
 %   func_3rd_deriv = @(x) x.^2 .* (3*fitCoef3(1)) + x .* (2*fitCoef3(2)) + fitCoef3(3);
%   func_length = @(x) sqrt( abs( 3 .* x.^2 .* fitCoef3(1) + 2 .* x .* fitCoef3(2) + fitCoef3(3) + 1 ) );

    % [y,x] are pixels in Lc image where this cell is located
    [y,x] = find(Lc == cellnum); % note: returns (row, column), which will be used as (y,x)
    x_rot = x*cos(phi) - y*sin(phi); % mathematical rotation
    y_rot = x*sin(phi) + y*cos(phi); % mathematical rotation
    
    xpixlength=max(x_rot)-min(x_rot);
    
    if xpixlength <250; 
         % function f(x) describing bacterial mid-axis
            func_3rd = @(x)                  x.^3 .* fitCoef3(1) +      x.^2 .* fitCoef3(2) + x .* fitCoef3(3)  + fitCoef3(4);
            
         % function l(x) = sqrt(1+f'(x)^2), to be integrated to determine length   
         func_length = @(x) sqrt( (abs( 3 .* x.^2 .* fitCoef3(1) + 2 .* x .* fitCoef3(2)    + fitCoef3(3)).^2   + 1 ) ); 
         
        order=3;
    else
        % 7TH DEGREE FIT
        
        % function f(x) describing bacterial mid-axis
           func_3rd = @(x)                  x.^7 .* fitCoef7(1) +     x.^6 .* fitCoef7(2) +     x.^5 .* fitCoef7(3) +     x.^4 .* fitCoef7(4)       +     x.^3 .* fitCoef7(5) +    x.^2 .* fitCoef7(6) + x .* fitCoef7(7) + fitCoef7(8);% RUTGER ADDED
        % function l(x) = sqrt(1+f'(x)^2), to be integrated to determine length   
        func_length = @(x) sqrt( (abs( 7 .* x.^6 .* fitCoef7(1) + 6 .*x.^5 .* fitCoef7(2) + 5 .*x.^4 .* fitCoef7(3) + 4 .*x.^3  * fitCoef7(4)       + 3 .*x.^2  * fitCoef7(5) +2 .*x.^1  * fitCoef7(6) + fitCoef7(7)).^2 + 1) ); % RUTGER ADDED
        
        order=7;
    end 

%  if s==1000;
%     s
% end
    schnitzcells(s).x_rot= x_rot;
    schnitzcells(s).y_rot= y_rot;
    schnitzcells(s).order= order;
%TRIES TO OBTAIN RIGHT CURVATURE GRAPH
%      y_rot_new1=func_3rd(x_rot);
%      figure;plot(x_rot,y_rot_new1);
%     figure;plot(x_rot,y_rot_new1);
%     curv=quick_curv(x_rot,y_rot_new1);
%     figure;plot(x_rot_intersect,y_rot_intersect);
%     figure;plot(x_rot,y_line_rot_left);
%     figure;plot(x_rot_line_left,y_line_rot_left);



    %--------------------------------------------------------------------
    % Determine length from 3rd degree polynomial, IMPROVED OLD WAY 
    % Adding: fitCoef3b_x_rot_left, fitCoef3b_x_rot_right, length_fitCoef3b
    %--------------------------------------------------------------------    
    % determine for intermediate x_rot's, what y_rot the fitted line has
    y_line_rot_left = func_3rd(x_rot-0.5);
    y_line_rot_right = func_3rd(x_rot+0.5);

    % make sure that y_line_rot_right contains lowest value
    for i = 1:length(y_line_rot_left)
      if y_line_rot_left(i) < y_line_rot_right(i)
        temp = y_line_rot_left(i);
        y_line_rot_left(i) = y_line_rot_right(i);
        y_line_rot_right(i) = temp;
      end
    end
    
    % Check which y_rot is closer than 0.5 away from line
    idx_pixel_intersects = find( y_rot >= y_line_rot_right-0.5 & y_rot <= y_line_rot_left+0.5 );
    x_rot_intersect = x_rot(idx_pixel_intersects);
    y_rot_intersect = y_rot(idx_pixel_intersects);

    % Determine length
    x_rot_intersect_left  = min(x_rot_intersect);
    x_rot_intersect_right = max(x_rot_intersect);
    schnitzcells(s).fitCoef3b_x_rot_left(age)   = x_rot_intersect_left;
    schnitzcells(s).fitCoef3b_x_rot_right(age)  = x_rot_intersect_right;
    schnitzcells(s).length_fitCoef3b(age) = quad(func_length, x_rot_intersect_left,x_rot_intersect_right) * p.micronsPerPixel;


   %zz=medfilt1(y_rot_intersect,21);zz1=smooth(zz,2);figure;plot(x_rot_intersect,zz1);
   % ynew_length
    %xnew_length
    %--------------------------------------------------------------------
    % Determine length from 3rd degree polynomial, NEW WAY
    % Adding: fitNew_x_rot_left, fitNew_x_rot_right, length_fitNew
    %--------------------------------------------------------------------    
    x_rot_line_left   = min(x_rot)-2:0.1:min(x_rot)+6; 
    x_rot_line_right  = max(x_rot)-6:0.1:max(x_rot)+2; 
    y_line_rot_left   = func_3rd(x_rot_line_left);
    y_line_rot_right  = func_3rd(x_rot_line_right);
 
    total_dist_25_points_left   = []; 
    total_dist_25_points_right  = []; 
    left_set = 0; right_set = 0;
    for i = 1:length(x_rot_line_left)
      dist_squared_left = (x_rot-x_rot_line_left(i)).^2 + (y_rot-y_line_rot_left(i)).^2;
      dist_squared_left_sorted = sort(dist_squared_left);
      total_dist_25_points_left(i) = sum(dist_squared_left_sorted(1:25)) * 0.01;
      
      dist_squared_right = (x_rot-x_rot_line_right(length(x_rot_line_left)+1-i)).^2 + (y_rot-y_line_rot_right(length(x_rot_line_left)+1-i)).^2;
      dist_squared_right_sorted = sort(dist_squared_right);
      total_dist_25_points_right(i) = sum(dist_squared_right_sorted(1:25)) * 0.01;
      
      if total_dist_25_points_left(i) < 1.1 & ~left_set
        schnitzcells(s).fitNew_x_rot_left(age)   = x_rot_line_left(i);
        left_set = 1;
        if i==1
          disp(['error: left: in fr' str3(currFrameNum) ' cellno ' str3(cellnum) ' schnitz ' str3(s) ]);
        end
      end
      if total_dist_25_points_right(i) < 1.1 & ~right_set
        schnitzcells(s).fitNew_x_rot_right(age)   = x_rot_line_right(length(x_rot_line_left)+1-i);
        right_set = 1;
        if i==1
          disp(['error: right: in fr' str3(currFrameNum) ' cellno ' str3(cellnum) ' schnitz ' str3(s) ]);
        end
      end
    end
    
    if ~left_set
      disp(['error: total_dist_25_points_left did not set left: in fr' str3(currFrameNum) ' cellno ' str3(cellnum) ' schnitz ' str3(s) ]);
      schnitzcells(s).fitNew_x_rot_left(age) = x_rot_line_left(end);
    end
    if ~right_set
      disp(['error: total_dist_25_points_left did not set right: in fr' str3(currFrameNum) ' cellno ' str3(cellnum) ' schnitz ' str3(s) ]);
      schnitzcells(s).fitNew_x_rot_right(age) = x_rot_line_right(1);
    end
    schnitzcells(s).length_fitNew(age) = quad(func_length, schnitzcells(s).fitNew_x_rot_left(age),schnitzcells(s).fitNew_x_rot_right(age)) * p.micronsPerPixel;
    
    %%% Find curvature
    
    x_values=min(schnitzcells(s).fitNew_x_rot_left(age)):0.5:max(schnitzcells(s).fitNew_x_rot_right(age));
    y_values=func_3rd(x_values);
    schnitzcells(s).x_values_end= x_values;
    schnitzcells(s).y_values_end= y_values;
    %figure(10); plot(x_values,y_values); pause; close figure 10
    %curv=quick_curv(x_values,y_values);
    %schnitzcells(s). curv= curv;
    s;
    framenum;

    if p.onScreen ~= 2 & find(p.schnitzNum==s) % set range to output as fig with p.schnitzNum
      %--------------------------------------------------------------------
      % calculate some data
      %--------------------------------------------------------------------
      size_Lc = size(Lc);
      
      % calc extended fitted line
      x_rot_line_total   = min(x_rot)-12:max(x_rot)+12; 
      y_rot_line_total   = func_3rd(x_rot_line_total);

      % Rotate extended fitted line back and make pixels
      x_line_total = x_rot_line_total*cos(-phi) - y_rot_line_total*sin(-phi); % mathematical rotation
      y_line_total = x_rot_line_total*sin(-phi) + y_rot_line_total*cos(-phi); % mathematical rotation

      % make sure that x and y are within Lc
      idx = find( round(x_line_total) <= size_Lc(2) & round(x_line_total) > 0 & round(y_line_total) <= size_Lc(1) & round(y_line_total) > 0 );
      x_line_total_good = x_line_total( idx );
      y_line_total_good = y_line_total( idx );
      
      line_total_idx = sub2ind(size_Lc, round(y_line_total_good), round(x_line_total_good)) ;
      line_total = zeros(size_Lc);
      line_total( line_total_idx ) = 1;

      % Calc distance to seg points, and phase contrast alonf extended fitted line
      total_dist_25_points_total    = [];
      phaseContrast_total = [];
      for i = 1:length(x_rot_line_total)
        dist_squared_total = (x_rot-x_rot_line_total(i)).^2 + (y_rot-y_rot_line_total(i)).^2;
        dist_squared_total_sorted = sort(dist_squared_total);
        total_dist_25_points_total(i) = sum(dist_squared_total_sorted(1:25)) * 0.01;
        
        phaseContrast_total(i) = interp2(phsub,x_line_total(i),y_line_total(i));
      end

      % normalize phaseContrast
      norm_idx = [round(length(phaseContrast_total)/3):round(2*length(phaseContrast_total)/3)];
      phaseContrast_total = (0.3/median(phaseContrast_total(norm_idx))) * phaseContrast_total;

      % make seg image
      Lc_cell = zeros(size(Lc));
      Lc_cell(find(Lc==cellnum)) = 1;
      Lc_cell(find(Lc~=cellnum & Lc~=0)) = 2;

      % get sub part of phase contrast
      cell = zeros(size(Lc));
      cell(find(Lc==cellnum)) = 1;
      cell = +(cell > 0);
      [fx,fy] = find(cell);
      xmin = max(min(fx)-5,1);
      xmax = min(max(fx)+5,size(cell,1));
      ymin = max(min(fy)-5,1);
      ymax = min(max(fy)+5,size(cell,2));
      subphcell = phsub(xmin:xmax, ymin:ymax);
      
      % autocontrast of phase
      [subphcell_index, map_gray] = gray2ind(imadjust(subphcell), 64);

      % add fitted line to phase contrast
      map_size = size(map_gray);
      map_gray(map_size(1)+1,:) = [ 0 0 1 ];
      subline_total = line_total(xmin:xmax, ymin:ymax);
      subphcell_index(find(subline_total>0)) = map_size(1);
      
      %--------------------------------------------------------------------
      % Make figure
      %--------------------------------------------------------------------    
      scrsz = get(0, 'ScreenSize');
      figureName = [p.movieDate ' ' p.movieName ' schnitz ' str4(s) ' frame ' str3(schnitzcells(s).frame_nrs(age))]; % MW 2014/06/11 removal N+1 bug
      figureFileName = ['length_schnitz' str4(s) 'frame' str3(schnitzcells(s).frame_nrs(age))]; % MW 2014/06/11 removal N+1 bug
      fig1 = figure('Position', [151 scrsz(4)-200 scrsz(3)-130 scrsz(4)-200], 'Name', figureName, 'visible','off'); %[1 scrsz(4)-200 scrsz(3) scrsz(4)-200]

      %--------------------------------------------------------------------
      % plot phase with fitted line
      %--------------------------------------------------------------------    
      subplot(2,3,1);
      subimage(subphcell_index, map_gray);
      axis off;
      
      % add title
      title(figureName,'interpreter','none','FontWeight','bold','FontSize',10);
      
      %--------------------------------------------------------------------
      % plot segmentation of cell
      %--------------------------------------------------------------------    
      subplot(2,3,4);
      DJK_imshowlabel(Lc_cell,'randomize',0);
      
      %--------------------------------------------------------------------
      % plot rotated x & y
      %--------------------------------------------------------------------    
      subplot(2,3,[2 3]);
      xlim([min(x_rot_line_total) max(x_rot_line_total)]);
      set(gca,'xtick',[],'ytick',[])

      % plot length OLD WAY
      hold on; plot([schnitzcells(s).fitCoef3_x_rot_left(age) schnitzcells(s).fitCoef3_x_rot_left(age)], [min(-y_rot) max(-y_rot)], 'g-', 'LineWidth',1);
      hold on; plot([schnitzcells(s).fitCoef3_x_rot_right(age) schnitzcells(s).fitCoef3_x_rot_right(age)], [min(-y_rot) max(-y_rot)], 'g-', 'LineWidth',1);

      % plot length OLD WAY b
      hold on; plot([schnitzcells(s).fitCoef3b_x_rot_left(age) schnitzcells(s).fitCoef3b_x_rot_left(age)], [min(-y_rot)+2 max(-y_rot)-2], 'g-', 'LineWidth',3);
      hold on; plot([schnitzcells(s).fitCoef3b_x_rot_right(age) schnitzcells(s).fitCoef3b_x_rot_right(age)], [min(-y_rot)+2 max(-y_rot)-2], 'g-', 'LineWidth',3);

      % plot length NEW WAY
      hold on; plot([schnitzcells(s).fitNew_x_rot_left(age) schnitzcells(s).fitNew_x_rot_left(age)], [min(-y_rot)+3 max(-y_rot)-3], 'b-', 'LineWidth',3);
      hold on; plot([schnitzcells(s).fitNew_x_rot_right(age) schnitzcells(s).fitNew_x_rot_right(age)], [min(-y_rot)+3 max(-y_rot)-3], 'b-', 'LineWidth',3);

      % plot fitted line
      hold on; plot(x_rot_line_total,-y_rot_line_total, 'b-', 'LineWidth',5);

      % plot rotated x & y
      hold on; plot(x_rot,-y_rot, 'k.');
    
      %--------------------------------------------------------------------
      % plot calculated length data
      %--------------------------------------------------------------------    
      subplot(2,3,[5 6]);
      xlim([min(x_rot_line_total) max(x_rot_line_total)]);
      ylim([0 5]);

      % plot length OLD WAY
      hold on; plot([schnitzcells(s).fitCoef3_x_rot_left(age) schnitzcells(s).fitCoef3_x_rot_left(age)], [0.1 1.9], 'g-', 'LineWidth',1);
      hold on; plot([schnitzcells(s).fitCoef3_x_rot_right(age) schnitzcells(s).fitCoef3_x_rot_right(age)], [0.1 1.9], 'g-', 'LineWidth',1);

      % plot length OLD WAY b
      hold on; plot([schnitzcells(s).fitCoef3b_x_rot_left(age) schnitzcells(s).fitCoef3b_x_rot_left(age)], [0.3 1.7], 'g-', 'LineWidth',3);
      hold on; plot([schnitzcells(s).fitCoef3b_x_rot_right(age) schnitzcells(s).fitCoef3b_x_rot_right(age)], [0.3 1.7], 'g-', 'LineWidth',3);

      % plot length NEW WAY
      hold on; plot([schnitzcells(s).fitNew_x_rot_left(age) schnitzcells(s).fitNew_x_rot_left(age)], [0.5 1.5], 'b-', 'LineWidth',3);
      hold on; plot([schnitzcells(s).fitNew_x_rot_right(age) schnitzcells(s).fitNew_x_rot_right(age)], [0.5 1.5], 'b-', 'LineWidth',3);

      % plot distance seg
      hold on; plot(x_rot_line_total,total_dist_25_points_total, 'k-', 'LineWidth',1);

      % plot phase
      hold on; plot(x_rot_line_total,phaseContrast_total, 'r-', 'LineWidth',1);

      % calc derivative of phase and plot *5 + 1
      d_phaseContrast_total = []; d_x_rot_line_total = [];
      for i = 1:length(phaseContrast_total)-1
        d_x_rot_line_total(i) = x_rot_line_total(i) + 0.5*(x_rot_line_total(i+1) - x_rot_line_total(i));
        d_phaseContrast_total(i) = 5 * (phaseContrast_total(i+1) - phaseContrast_total(i)) + 1;
      end
      hold on; plot(d_x_rot_line_total,d_phaseContrast_total, 'r:', 'LineWidth',1);
      
      %--------------------------------------------------------------------------
      % ASKING TO SAVE FIGURE
      %--------------------------------------------------------------------------
      % Ask to save the figure
      if p.onScreen
        set(fig1,'visible','on');
        saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','Yes and Close','No','Yes');
        pause(0.2);
      else
        saveFigInput='Yes and Close';
      end
      if (upper(saveFigInput(1))=='Y')
          saveSameSize(fig1,'file',[p.DJK_saveDir figureFileName '.png'], 'format', 'png');
          if (strcmp(saveFigInput,'Yes and Close'))
              close(fig1);
              pause(0.2);
          end
        disp([' * Saved plot in ' figureFileName '.png']);
      end
      %--------------------------------------------------------------------------    
    end
  end % loop over schnitzesForFrame
end % loop over frames
%--------------------------------------------------------------------------



for k=1:max(size(schnitzcells));
 x_values=schnitzcells(k).x_values_end(1,:);
 y_values=schnitzcells(k).y_values_end(1,:);
 curv_all=quick_curv(x_values,y_values);
 schnitzcells(k). curv= curv_all;
end
%--------------------------------------------------------------------------
% Save extacted data
%--------------------------------------------------------------------------
save(p.schnitzName, 'schnitzcells');
disp(['Save in ''' p.schnitzName ''' completed...']);
%--------------------------------------------------------------------------
