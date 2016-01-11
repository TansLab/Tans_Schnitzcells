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

% Last Modified by GUIDE v2.5 08-Jan-2016 18:46:02

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
% uiwait(handles.figure1);


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
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections

% call script, runsections defines which sectino of sript to run
runsections = 'reloadfile'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'createpfull'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'cropimagesfull'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'segmentation'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'manualchecksegfull'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'quickanalysis'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'trackandmanualcorrections'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
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

% update workspace

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
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'createbackup'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'correctionsandanalysis'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)
assignin ('base','schnitzcells',schnitzcells)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'loadpforcropped'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global runsections 

% call script, runsections defines which sectino of sript to run
runsections = 'makeoutputfull'
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings)
assignin ('base','p',p)
assignin ('base','schnitzcells',schnitzcells)
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

global runsections customFrameRange p settings

% retrieve settings var
settings = evalin('base', 'settings');
p.overwrite=1;

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(settings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
customFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
settings.specialtracker = 'MW';
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 
% reset specialtracker
settings = rmfield(settings,'specialtracker');

% update workspace
assignin ('base','settings',settings);
assignin ('base','p',p);

% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global runsections customFrameRange p settings

% retrieve settings var
settings = evalin('base', 'settings');

% Set custom tracker:
settings.specialtracker = 'NW';
p.overwrite=1;

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(settings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
customFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
settings.specialtracker = 'NW';
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 
% reset specialtracker
settings = rmfield(settings,'specialtracker');

% update workspace
assignin ('base','settings',settings);
assignin ('base','p',p);

% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global runsections customFrameRange p settings

% retrieve settings var
settings = evalin('base', 'settings');
p.overwrite=1;

% Get custom framerange
prompt={'Enter frame range:'};
name = 'Range:';
defaultans = {mat2str(settings.currentFrameRange)};
answer = inputdlg(prompt,name,1,defaultans);
customFrameRange = str2num(answer{1});

% call script, runsections defines which sectino of sript to run
runsections = 'customtrackersoncustomrange'
% run w. default tracker
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings);
assignin ('base','p',p);


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Re-run segmentation for certain frame

global runsections customFrameRange p settings

% retrieve settings var
settings = evalin('base', 'settings');
p.overwrite=1;

% Get custom parameters
% Formulate questions
prompt={'Frame number or range:','Slices','Laplace of Gaussian smoothing:','Minimum depth:','Minimum cell area'};
defaultanswers = {'1','[1,2,3]','2','5','250'}
name = 'Enter segmentation parameters';
% Get answers using prompt
answers = inputdlg(prompt,name,1,defaultanswers);
% Rename answers and pass through using settings
settings.TOREDOFRAME            = str2num(answers{1});
settings.SLICESTEMPORARY        = str2num(answers{2});
settings.LOGSMOOTHINGTEMPORARY  = str2num(answers{3});
settings.MINDEPTHTEMP           = str2num(answers{4});
settings.MINCELLAREA            = str2num(answers{5});

% call script, runsections defines which sectino of sript to run
runsections = 'redosegforframe'
% run w. default tracker
MW_analysis_attempt2_matlabinsteadexcel_plusGUI 

% update workspace
assignin ('base','settings',settings);
assignin ('base','p',p);
