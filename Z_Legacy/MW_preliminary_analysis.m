
% Script to perform preliminary analysis. 
% MW 2015/07
% 
% Replaces excel file


%% USER PARAMETERS
% ===

USEFULLIMAGE        = 1;
OVERWRITE           = 0;

DATE                = '2015-06-12';
SETUP               = 'setup1';
SOFTWAREPACKAGE     = 'metamorph';
CAMERA              = 'hamamatsu';
POSITIONNAME        = 'pos4';
FOLDERNAME          = 'T:\TRAVELING_DATA';
FRAMERANGE          = [16:5:246];
FLUOR1              = 'g';
FLUOR2              = 'none';
FLUOR3              = 'none';
	
LEFTTOP             = [373 5];
RIGHTBOTTOM         =	[2042 1768];
	
SLICES              = [1 2 3];
RANGEFILTSIZE       =	35;
MASKMARGIN          = 20;
LOGSMOOTHING        =	2;
MINCELLAREA         = 250;
GAUSSIANFILTER      = 5;
MINDEPTH            = 5;
NECKDEPTH           = 2;
	
XLIMIT              = [0 1000]; % for time axes
	
IMAGESPERMU         = 15;
SCHNITZESTOBEREMOVED = [];
AUTOFLUORESCENCE    = 0.94;

% Some additional parameters
STRAINNAME          = 'e.coli.amolf';

disp('Parameters set.');

%% Cropping

% determine crop area
p1 = DJK_initschnitz(POSITIONNAME,DATE,STRAINNAME,'rootDir',FOLDERNAME, 'cropLeftTop',LEFTTOP, 'cropRightBottom',RIGHTBOTTOM,'fluor1',FLUOR1,'fluor2',FLUOR2,'fluor3',FLUOR3,'setup',SETUP,'softwarePackage',SOFTWAREPACKAGE,'camera',CAMERA);
p1.imageDir=[FOLDERNAME '\' DATE '\' POSITIONNAME '\']
[LEFTTOP, RIGHTBOTTOM] = MW_determinecroparea(p1, FRAMERANGE)

% set up p1
p1 = DJK_initschnitz(POSITIONNAME,DATE,STRAINNAME,'rootDir',FOLDERNAME, 'cropLeftTop',LEFTTOP, 'cropRightBottom',RIGHTBOTTOM,'fluor1',FLUOR1,'fluor2',FLUOR2,'fluor3',FLUOR3,'setup',SETUP,'softwarePackage',SOFTWAREPACKAGE,'camera',CAMERA);
p1.imageDir=[FOLDERNAME '\' DATE '\' POSITIONNAME '\']

% ask user to continue
aString = input('Continue (y/n enter) ? >>', 's')
if ~(aString=='y')
    error('Analysis cancelled.');
end

disp('Continuing');

DJK_cropImages_3colors(p1, FRAMERANGE, LEFTTOP, RIGHTBOTTOM, 'cropName', 'pos4crop');

p1 = DJK_initschnitz('pos4crop','2015-06-12','e.coli.amolf','rootDir','T:\TRAVELING_DATA\', 'cropLeftTop',[373 5], 'cropRightBottom',[2042 1768],'fluor1','g','fluor2','none','fluor3','none','setup','setup1','softwarePackage','metamorph','camera','hamamatsu');

p1.useFullImage=USEFULLIMAGE;
p1.overwrite=OVERWRITE;


%%



