function varargout = gmtOlsa_gui_testClosed(varargin)
% GMTOLSA_GUI_TESTCLOSED MATLAB code for gmtOlsa_gui_testClosed.fig
%      GMTOLSA_GUI_TESTCLOSED, by itself, creates a new GMTOLSA_GUI_TESTCLOSED or raises the existing
%      singleton*.
%
%      H = GMTOLSA_GUI_TESTCLOSED returns the handle to a new GMTOLSA_GUI_TESTCLOSED or the handle to
%      the existing singleton*.
%
%      GMTOLSA_GUI_TESTCLOSED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GMTOLSA_GUI_TESTCLOSED.M with the given input arguments.
%
%      GMTOLSA_GUI_TESTCLOSED('Property','Value',...) creates a new GMTOLSA_GUI_TESTCLOSED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gmtOlsa_gui_testClosed_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gmtOlsa_gui_testClosed_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gmtOlsa_gui_testClosed

% Last Modified by GUIDE v2.5 23-Aug-2021 14:19:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gmtOlsa_gui_testClosed_OpeningFcn, ...
                   'gui_OutputFcn',  @gmtOlsa_gui_testClosed_OutputFcn, ...
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


% --- Executes just before gmtOlsa_gui_testClosed is made visible.
function gmtOlsa_gui_testClosed_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gmtOlsa_gui_testClosed (see VARARGIN)

% Choose default command line output for gmtOlsa_gui_testClosed
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gmtOlsa_gui_testClosed wait for user response (see UIRESUME)
% uiwait(handles.User_interact_mask);


% --- Outputs from this function are returned to the command line.
function varargout = gmtOlsa_gui_testClosed_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
    set(handles.bg_control, 'Userdata', get(hObject, 'Userdata'));
    uiresume;
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pb_cancel


% --- Executes during object creation, after setting all properties.
function User_interact_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to User_interact_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in bg_control.
function bg_control_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bg_control 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
    


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over togglebutton5.
function togglebutton5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on togglebutton5 and none of its controls.
function togglebutton5_KeyPressFcn(hObject, eventdata, handles)
   
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in bg_col1.
function bg_col1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bg_col1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pb_next and none of its controls.
function pb_next_KeyPressFcn(hObject, eventdata, handles)   

% hObject    handle to pb_next (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pb_cancel and none of its controls.
function pb_cancel_KeyPressFcn(hObject, eventdata, handles)    
% hObject    handle to pb_cancel (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pb_replay and none of its controls.
function pb_replay_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pb_replay (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_next.
function pb_next_Callback(hObject, eventdata, handles)
    set(handles.bg_control, 'UserData', get(hObject, 'UserData'));
    uiresume;
% hObject    handle to pb_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_replay.
function pb_replay_Callback(hObject, eventdata, handles)
    set(handles.bg_control, 'UserData', get(hObject, 'UserData'));
    uiresume;
% hObject    handle to pb_replay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
