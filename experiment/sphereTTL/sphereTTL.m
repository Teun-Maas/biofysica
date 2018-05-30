function varargout = sphereTTL(varargin)
% SPHERETTL MATLAB code for sphereTTL.fig
%      SPHERETTL, by itself, creates a new SPHERETTL or raises the existing
%      singleton*.
%
%      H = SPHERETTL returns the handle to a new SPHERETTL or the handle to
%      the existing singleton*.
%
%      SPHERETTL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPHERETTL.M with the given input arguments.
%
%      SPHERETTL('Property','Value',...) creates a new SPHERETTL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sphereTTL_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sphereTTL_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sphereTTL

% Last Modified by GUIDE v2.5 30-May-2018 14:38:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @sphereTTL_OpeningFcn, ...
	'gui_OutputFcn',  @sphereTTL_OutputFcn, ...
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


% --- Executes just before sphereTTL is made visible.
function sphereTTL_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sphereTTL (see VARARGIN)


%% Get all important trial and configuration info
handles			= getcfgdefaults(handles);
sethandles(handles);
handles			= gethandles(handles);

%% Initialize TDT
handles			= tdt_init(handles); % Turning on, Loading, and Running Circuits

%% Start experiment
handles.data	= struct([]);
handles			= setupShow(handles);

% Choose default command line output for sphereTTL
handles.output	= hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sphereTTL wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = sphereTTL_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in checkbox_RA16connect.
function checkbox_RA16connect_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to checkbox_RA16connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RA16connect

% --- Executes on button press in checkbox_RA16load.
function checkbox_RA16load_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RA16load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RA16load

% --- Executes on button press in checkbox_RA16run.
function checkbox_RA16run_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RA16run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RA16run

% --- Executes on button press in checkbox_RP21connect.
function checkbox_RP21connect_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP21connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP21connect


% --- Executes on button press in checkbox_RP21load.
function checkbox_RP21load_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP21load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP21load


% --- Executes on button press in checkbox_RP21run.
function checkbox_RP21run_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP21run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP21run


% --- Executes on button press in checkbox_RP22connect.
function checkbox_RP22connect_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP22connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP22connect

% --- Executes on button press in checkbox_RP22load.
function checkbox_RP22load_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP22load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP22load

% --- Executes on button press in checkbox_RP22run.
function checkbox_RP22run_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RP22run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RP22run

% --------------------------------------------------------------------
function experiment_Callback(hObject, eventdata, handles)
% hObject    handle to experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
runExperiment(handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
runExperimentTTL(handles)

% --------------------------------------------------------------------
function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function edit_subjectid_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subjectid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subjectid as text
%        str2double(get(hObject,'String')) returns contents of edit_subjectid as a double
handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_subjectid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subjectid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_expinitials.
function popupmenu_expinitials_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_expinitials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_expinitials contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_expinitials
handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_expinitials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_expinitials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function edit_date_Callback(hObject, eventdata, handles)
% hObject    handle to edit_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_date as text
%        str2double(get(hObject,'String')) returns contents of edit_date as a double
handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function edit_block_Callback(hObject, eventdata, handles)
% hObject    handle to edit_block (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_block as text
%        str2double(get(hObject,'String')) returns contents of edit_block as a double
handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_block_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_block (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_exp.
function popupmenu_exp_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_exp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_exp
handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_exp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_cfg.
function popupmenu_cfg_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cfg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cfg

handles = gethandles(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_cfg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end
