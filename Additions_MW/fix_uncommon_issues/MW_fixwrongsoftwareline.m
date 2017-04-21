

% This fixes a very specific problem.
%
% In dataset below, one image was opened and saved with another software
% program, messing up schnitzcell's ability to recognize the correct
% software.
%
% Solution: let matlab open and save the image, adding the description with
% the info manually.

myDescription = ['Exposure: 150 ms' 10 'Binning: 2 x 2' 10 'Region: 1392 x 1040, offset at (0, 0)' 10 'Subtract: Off' 10 'Shading: Off' 10 'Digitizer: 10 MHz' 10 'Gain: Gain 2 (4x)' 10 'Camera Shutter: Always Open' 10 'Clear Count: 2' 10 'Clear Mode: CLEAR PRE SEQUENCE' 10 'Frames to Average: 1' 10 'Trigger Mode: Normal (TIMED)' 10 'Temperature: -5.05' 10 'ANNOTATION_DAAN' 10 '' 10 'Camera.Digital.Exposure = 150' 10 '' 10 'StagePosition.StageX = -1263.9' 10 'StagePosition.StageY = 3124' 10 'StagePosition.StageZ = 4503.8' 10 '' 10 'Device.Illumination.Setting = FL - engfp' 10 'Device.Magnification.Setting = PlanFluor 100X' 10 'Device.Stage.XPosition = -1263.9' 10 'Device.Stage.YPosition = 3124' 10 'Device.Focus.CurPos = 4503.8' 10 '' 10 'Image.StageX = -1263.9' 10 'Image.StageY = 3124' 10 'Image.ZAbsolute = 4503.8' 10 'Image.Time.Create.Date = 2.45656e+006' 10 'Image.Time.Create.time = 2.07302e+007' 10 'Image.IllumSetting = FL - mCherry' 10 'Image.MagSetting = PlanFluor 100X' 10 'Image.NA = 1.4' 10 'Image.Name = UntitledDateTime: 2013:09:25 05:45:34Software: MetaMorph 7.1.4.0'];

myImg=imread('G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\images\pos4crop-g-1244_original.tif');
imwrite(myImg,'G:\EXPERIMENTAL_DATA_2016\a_incoming\2013-09-24_Rutger1_TET\pos4crop\images\pos4crop-g-1244_fixed.tif','Description',myDescription);

