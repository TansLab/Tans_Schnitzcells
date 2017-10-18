function varargout = MW_GUI_schnitzcells(varargin)
% MW_GUI_SCHNITZCELLS MATLAB code for MW_GUI_schnitzcells.fig
%
% Note: type 
% >> guide MW_GUI_schnitzcells 
% To easily edit.
%
%      MW_GUI_SCHNITZCELLS, by itself, creates a new MW_GUI_SCHNITZCELLS or raises the existing
%      singleton*.
%
%      H = MW_GUI_SCHNITZCELLS returns the handle to a new MW_GUI_SCHNITZCELLS or the handle to
%      the existing singleton*.
%
%      MW_GUI_SCHNITZCELLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MW_GUI_SCHNITZCELLS.M with the given input arguments.
%
%      MW_GUI_SCHNITZCELLS('Property','Value',...) creates a new MW_GUI_SCHNITZCELLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MW_GUI_schnitzcells_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MW_GUI_schnitzcells_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MW_GUI_schnitzcells

% Last Modified by GUIDE v2.5 01-Jun-2017 17:15:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MW_GUI_schnitzcells_OpeningFcn, ...
                   'gui_OutputFcn',  @MW_GUI_schnitzcells_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MW_GUI_schnitzcells is made visible.
function MW_GUI_schnitzcells_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MW_GUI_schnitzcells (see VARARGIN)

% Choose default command line output for MW_GUI_schnitzcells
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MW_GUI_schnitzcells wait for user response (see UIRESUME)
% uiwait(handles.schnitzcellsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = MW_GUI_schnitzcells_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections

% call script, runsections defines which sectino of sript to run
runsections = 'loadfile'
Schnitzcells_masterscript

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','alldata',alldata)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections

% call script, runsections defines which sectino of sript to run
runsections = 'reloadfile'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','alldata',alldata)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'createp'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% retrieve alldata var from base
alldata = evalin('base', 'alldata');

% call script, runsections defines which sectino of sript to run
runsections = 'cropimages'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% continue as usual
global runsections p

% clean up workspace from crop step
evalin('base', 'clear alldata' );

% Make sure settings are obtained
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

% Pass on whether we want to use full image for segmentation (default = 0)
p.useFullImage = get(handles.checkbox4,'Value');

% call script, runsections defines which sectino of sript to run
runsections = 'segmentation'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections ourSettings p
% obtain ourSettings
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');
    
% call script, runsections defines which sectino of sript to run
runsections = 'manualchecksegfull'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections ourSettings

% call script, runsections defines which sectino of sript to run
runsections = 'quickanalysis'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'trackandmanualcorrections'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'analysispreliminary'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)
assignin ('base','schnitzcells',schnitzcells)

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run

% update workspace

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'trackpreliminary'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)



% call script, runsections defines which sectino of sript to run

% update workspace

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run

% update workspace

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run

% update workspace

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% run appropriate section
runsections = 'cropimagespreliminary'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)


% call script, runsections defines which sectino of sript to run

% update workspace

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run

% update workspace


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'makemovie'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'createbackup'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'correctionsandanalysis'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)
if exist('schnitzcells','var')
    assignin ('base','schnitzcells',schnitzcells)
else
    disp('(Schnitzcells not exported to workspace; it does get saved by subfunctions..)');
end


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'loadpforcropped'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'makeoutputfull'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)
assignin ('base','schnitzcells',schnitzcells)
assignin ('base','s_rm',s_rm)
assignin ('base','output',output)


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections p ourSettings

% retrieve ourSettings var
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');
p.overwrite = get(handles.checkbox5,'Value');

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(ourSettings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
ourSettings.retrackFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
ourSettings.specialtracker = 'MW';
Schnitzcells_masterscript 
% reset specialtracker
ourSettings = rmfield(ourSettings,'specialtracker');

% update workspace
assignin ('base','ourSettings',ourSettings);
assignin ('base','p',p);

% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global runsections customFrameRange p ourSettings

% retrieve ourSettings var
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

% Set custom tracker:
ourSettings.specialtracker = 'NW';
p.overwrite = get(handles.checkbox5,'Value');

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(ourSettings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
ourSettings.retrackFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
ourSettings.specialtracker = 'NW';
Schnitzcells_masterscript 
% reset specialtracker
ourSettings = rmfield(ourSettings,'specialtracker');

% update workspace
assignin ('base','ourSettings',ourSettings);
assignin ('base','p',p);

% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections p ourSettings

% retrieve ourSettings var
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');
p.overwrite = get(handles.checkbox5,'Value');

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(ourSettings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
ourSettings.retrackFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
% run w. default tracker
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings);
assignin ('base','p',p);


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Re-run segmentation for certain frame

global runsections p ourSettings

% retrieve our parameters
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');
p.overwrite=1;
% Pass on whether we want to use full image for segmentation (default = 0)
p.useFullImage = get(handles.checkbox4,'Value');

% Get custom parameters
% Formulate questions
prompt={'Frame number or range:','Slices','Laplace of Gaussian smoothing:','Minimum depth:','Minimum cell area','rangeFiltSize:', 'maskMargin:', 'GaussianFilter:', 'neckDepth:'};
defaultanswers = {'1','[1,2,3]','2','5','250','35','5','5','2'}
name = 'Enter segmentation parameters';
% Get answers using prompt
answers = inputdlg(prompt,name,1,defaultanswers);
% Rename answers and pass through using ourSettings
ourSettings.TOREDOFRAME            = str2num(answers{1});
ourSettings.SLICESTEMPORARY        = str2num(answers{2});
ourSettings.LOGSMOOTHINGTEMPORARY  = str2num(answers{3});
ourSettings.MINDEPTHTEMP           = str2num(answers{4});
ourSettings.MINCELLAREA            = str2num(answers{5});
ourSettings.RANGEFILTSIZE          = str2num(answers{6});
ourSettings.MASKMARGIN             = str2num(answers{7});
ourSettings.GAUSSIANFILTER         = str2num(answers{8});
ourSettings.NECKDEPTH              = str2num(answers{9});

% call script, runsections defines which sectino of sript to run
runsections = 'redosegforframe'
% run w. default tracker
Schnitzcells_masterscript 

p.overwrite=0;

% update workspace
assignin ('base','ourSettings',ourSettings);
assignin ('base','p',p);


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'checkaftercustom'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% retrieve ourSettings var from base
ourSettings = evalin('base', 'ourSettings');

% open config file
winopen(ourSettings.configfilepath);

% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ourSettings

% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');

% set flag
ourSettings.analysisType = 'preliminary';

% export flag
assignin ('base','ourSettings',ourSettings);
disp('Preliminary flag set.');

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ourSettings

% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');

% set flag
ourSettings.analysisType = 'full';

% export flag
assignin ('base','ourSettings',ourSettings)


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

myAnswer = questdlg(['Clear all parameters?'],'Confirmation required.','Yes','No','No');

if strcmp(myAnswer,'Yes');
   evalin('base', 'clear all;');
   disp('Workspace cleared.');
else
    disp('Wiping aborted.');
end



function slicenr_Callback(hObject, eventdata, handles)
% hObject    handle to slicenr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slicenr as text
%        str2double(get(hObject,'String')) returns contents of slicenr as a double


% --- Executes during object creation, after setting all properties.
function slicenr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicenr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function framenr_Callback(hObject, eventdata, handles)
% hObject    handle to framenr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framenr as text
%        str2double(get(hObject,'String')) returns contents of framenr as a double


% --- Executes during object creation, after setting all properties.
function framenr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framenr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

% input from user
framenr = get(handles.framenr,'String');
slicenr = get(handles.slicenr,'String');

myimg = imread([p.imageDir p.movieName '-p-' slicenr '-' framenr '.tif']);
h2=figure(); imshow(myimg, []);

h3=figure(); 
for i=1:3
    subplot(3,1,i); hold on;
    myimg = imread([p.imageDir p.movieName '-p-' num2str(i) '-' framenr '.tif']);
    imshow(myimg, []);
end


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

% input from user
framenr = get(handles.framenr,'String');
%slicenr = get(handles.slicenr,'String');

% crete fig and make (almost) full screen
h2=figure(); 
set(h2, 'units','normalized', 'Position', [0.05,0.05,.9,.8]); % left bottom width height

% plot all slices
for i=1:3
    subplottight(1,3,i); hold on;
    myimg = imread([p.imageDir p.movieName '-p-' num2str(i) '-' framenr '.tif']);
    imshow(myimg, []);
end



function framenrtracking_Callback(hObject, eventdata, handles)
% hObject    handle to framenrtracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framenrtracking as text
%        str2double(get(hObject,'String')) returns contents of framenrtracking as a double


% --- Executes during object creation, after setting all properties.
function framenrtracking_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framenrtracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

% input from user
framenr = get(handles.framenrtracking,'String');
%slicenr = get(handles.slicenr,'String');

MW_helperforlinkingframes(p, str2num(framenr));


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes during object creation, after setting all properties.
function pushbutton35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function framenr2_Callback(hObject, eventdata, handles)
% hObject    handle to framenr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framenr2 as text
%        str2double(get(hObject,'String')) returns contents of framenr2 as a double


% --- Executes during object creation, after setting all properties.
function framenr2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framenr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cellLabelMYNEEDLE_Callback(hObject, eventdata, handles)
% hObject    handle to cellLabelMYNEEDLE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cellLabelMYNEEDLE as text
%        str2double(get(hObject,'String')) returns contents of cellLabelMYNEEDLE as a double


% --- Executes during object creation, after setting all properties.
function cellLabelMYNEEDLE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellLabelMYNEEDLE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% retrieve ourSettings
ourSettings = evalin('base', 'ourSettings');
p = evalin('base', 'p');

%% input from user
framenr = get(handles.framenr2,'String');
cellLabelMYNEEDLE = str2num(get(handles.cellLabelMYNEEDLE,'String'));

%% load segfile
load([p.segmentationDir p.movieName 'seg' num2str(framenr) '.mat']);

%% make figure with highlighted cell
Lselect=Lc; Lselect(Lc>0)=1; Lselect(Lc==cellLabelMYNEEDLE)=2;
figure; clf; imshow(Lselect,[]);





% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global p

% retrieve ourSettings
p = evalin('base', 'p');
ourSettings = evalin('base', 'ourSettings');

myimg = imread([p.imageDir p.movieName '-p-' '2' '-' sprintf('%03d',min(ourSettings.currentFrameRange)) '.tif']);
h2=figure(); clf; imshow(myimg, []);

p.customColonyCenter = round(ginput(1))

close(h2);

%cellno1 = (location(2),location(1));

% update workspace
assignin ('base','p',p)



% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global p

% retrieve ourSettings
p = evalin('base', 'p');

winopen(p.dateDir);



% --- Executes on button press in skipFramesWithProblems.
function skipFramesWithProblems_Callback(hObject, eventdata, handles)
% hObject    handle to skipFramesWithProblems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skipFramesWithProblems

% set skipFramesWithProblems

% --- Executes on button press in showPhaseImage.
function showPhaseImage_Callback(hObject, eventdata, handles)
% hObject    handle to showPhaseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPhaseImage
% set showPhaseImage

% --- Executes on button press in showLargeCentroids.
function showLargeCentroids_Callback(hObject, eventdata, handles)
% hObject    handle to showLargeCentroids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showLargeCentroids

% set showLargeCentroids

% --- Executes on button press in dontShowExtendedReport.
function dontShowExtendedReport_Callback(hObject, eventdata, handles)
% hObject    handle to dontShowExtendedReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dontShowExtendedReport

% set dontShowExtendedReport


function mothermachine_Callback(hObject, eventdata, handles)
% hObject    handle to mothermachine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mothermachine as text
%        str2double(get(hObject,'String')) returns contents of mothermachine as a double


% --- Executes during object creation, after setting all properties.
function mothermachine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mothermachine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton52.
function pushbutton52_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% retrieve our parameters
p = evalin('base', 'p');

p.showPhaseImage = ...
    get(handles.showPhaseImage,'Value');
p.showLargeCentroids = ...
    get(handles.showLargeCentroids,'Value');
p.dontShowExtendedReport = ...
    get(handles.dontShowExtendedReport,'Value');
p.skipFramesWithProblems = ...
    get(handles.skipFramesWithProblems,'Value');
p.mothermachine = ...
    str2num(get(handles.mothermachine,'String'));

%{
showPhaseImage
showLargeCentroids
dontShowExtendedReport
ignoreFailedChecksTracker
mothermachine
%}

% update workspace
assignin ('base','p',p);

disp('p structure updated w. additional options.');




% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'makemovieraw'
Schnitzcells_masterscript 

% update workspace
assignin ('base','ourSettings',ourSettings)
assignin ('base','p',p)
