function varargout = matrix_gui_main(varargin)
% matrix_GUI_MAIN MATLAB code for matrix_gui_main.fig
%      matrix_GUI_MAIN, by itself, creates a new matrix_GUI_MAIN or raises the existing
%      singleton*.
%
%      H = matrix_GUI_MAIN returns the handle to a new matrix_GUI_MAIN or the handle to
%      the existing singleton*.
%
%      matrix_GUI_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in matrix_GUI_MAIN.M with the given input arguments.
%
%      matrix_GUI_MAIN('Property','Value',...) creates a new matrix_GUI_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matrix_gui_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matrix_gui_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matrix_gui_main

% Last Modified by GUIDE v2.5 02-Sep-2021 12:23:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matrix_gui_main_OpeningFcn, ...
                   'gui_OutputFcn',  @matrix_gui_main_OutputFcn, ...
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


% --- Executes just before matrix_gui_main is made visible.
function matrix_gui_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matrix_gui_main (see VARARGIN)
% Choose default command line output for matrix_gui_main
handles.output = hObject;

filename = 'default.mat';
path = matrix_getPathSettings;
matrix_loadSettings(handles, fullfile(path.SETTINGS, filename));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes matrix_gui_main wait for user response (see UIRESUME)
% uiwait(handles.fig_main);


% --- Outputs from this function are returned to the command line.
function varargout = matrix_gui_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_start.
function pb_start_Callback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pb_start.
function pb_start_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pb_start and none of its controls.
function pb_start_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_quit.
function pb_quit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Do you really want to quit?','Quit','Yes','No','Yes');
if strcmp(button, 'Yes')
    close(handles.fig_main);
end



function edt_speechDelay_Callback(hObject, eventdata, handles)
% hObject    handle to edt_speechDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_speechDelay as text
%        str2double(get(hObject,'String')) returns contents of edt_speechDelay as a double


% --- Executes during object creation, after setting all properties.
function edt_speechDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_speechDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_noiseLvl_Callback(hObject, eventdata, handles)
% hObject    handle to edt_noiseLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_noiseLvl as text
%        str2double(get(hObject,'String')) returns contents of edt_noiseLvl as a double


% --- Executes during object creation, after setting all properties.
function edt_noiseLvl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_noiseLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_crossDur_Callback(hObject, eventdata, handles)
% hObject    handle to edt_crossDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_crossDur as text
%        str2double(get(hObject,'String')) returns contents of edt_crossDur as a double


% --- Executes during object creation, after setting all properties.
function edt_crossDur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_crossDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_speechLvl_Callback(hObject, eventdata, handles)
% hObject    handle to edt_speechLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_speechLvl as text
%        str2double(get(hObject,'String')) returns contents of edt_speechLvl as a double


% --- Executes during object creation, after setting all properties.
function edt_speechLvl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_speechLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_sensIn_Callback(hObject, eventdata, handles)
% hObject    handle to edt_sensIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_sensIn as text
%        str2double(get(hObject,'String')) returns contents of edt_sensIn as a double


% --- Executes during object creation, after setting all properties.
function edt_sensIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_sensIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_listNr.
function pop_listNr_Callback(hObject, eventdata, handles)
% hObject    handle to pop_listNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_listNr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_listNr


% --- Executes during object creation, after setting all properties.
function pop_listNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_listNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_stratFile.
function pb_stratFile_Callback(hObject, eventdata, handles)
% hObject    handle to pb_stratFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] =  uigetfile({'*.m';'*.mat'},'Strategy File Selector');
if (filename)
    pathname = strrep(pathname, pwd, ['.']);
    set(handles.txt_stratFile, 'String', [pathname filename]);
end


function edt_inputId_Callback(hObject, eventdata, handles)
% hObject    handle to edt_inputId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_inputId as text
%        str2double(get(hObject,'String')) returns contents of edt_inputId as a double


% --- Executes during object creation, after setting all properties.
function edt_inputId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_inputId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edt_subjectId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_subjectId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_sensOut_Callback(hObject, eventdata, handles)
% hObject    handle to edt_sensOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_sensOut as text
%        str2double(get(hObject,'String')) returns contents of edt_sensOut as a double


% --- Executes during object creation, after setting all properties.
function edt_sensOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_sensOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_resultsDir.
function pb_resultsDir_Callback(hObject, eventdata, handles)
% hObject    handle to pb_resultsDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = matrix_getPathSettings;
cd(path.RESULTS);
foldername = uigetdir(pwd, 'Choose root directory for test results');
cd(path.MATRIX)
if foldername
    %foldername = strrep(foldername, pwd, ['.']);
    set(handles.txt_resultsDir, 'String', foldername);
end

% --- Executes on button press in pb_noiseFile.
function pb_noiseFile_Callback(hObject, eventdata, handles)
% hObject    handle to pb_noiseFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = matrix_getPathSettings;
cd(path.SOUNDS_NL);
[filename, pathname] =  uigetfile({'*.wav'},'Noise File');
cd(path.MATRIX);
if filename
    pathname = strrep(pathname, pwd, ['.']);
    set(handles.txt_noiseFile, 'String', [pathname filename]);
end

% --- Executes on button press in cb_guess.
function cb_guess_Callback(hObject, eventdata, handles)
% hObject    handle to cb_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_guess



function edt_stratVarName_Callback(hObject, eventdata, handles)
% hObject    handle to edt_stratVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_stratVarName as text
%        str2double(get(hObject,'String')) returns contents of edt_stratVarName as a double


% --- Executes during object creation, after setting all properties.
function edt_stratVarName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_stratVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_stratView.
function pb_stratView_Callback(hObject, eventdata, handles)
% hObject    handle to pb_stratView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stratFile = get(handles.txt_stratFile, 'String');

if ~isempty(stratFile)
    matrix_loadSavedStrategy(stratFile, 1);
else
    warndlg('Please specify a strategy file and variable name');
end



% --- Executes on button press in pb_load.
function pb_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = matrix_getPathSettings;
cd(path.SETTINGS);
[filename, pathname] =  uigetfile({'*.mat'},'Load OLSA settings');
cd(path.MATRIX);
if filename
    pathname = strrep(pathname, pwd, ['.']);
    try
        matrix_loadSettings(handles, [pathname, filename]);
    catch ex
        disp('An error occured while trying to load the settings:');
        disp(ex.message);
        warndlg('An error occured while trying to load the settings.','Warning','modal');
    end
end



% --- Executes on button press in pb_save.
function pb_save_Callback(hObject, eventdata, handles)
% hObject    handle to pb_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = matrix_getPathSettings;
cd(path.SETTINGS);
[filename, pathname] =  uiputfile({'*.mat'},'Save OLSA settings');
cd(path.MATRIX);
if filename
    pathname = strrep(pathname, pwd, '.');
    try
        matrix_saveSettings(handles, [pathname, filename]);
    catch ex
        disp('An error occured while trying to save the settings:');
        disp(ex.message);
        warndlg('An error occured while trying to save the settings.','Warning','modal');
    end
end



function edt_outputId_Callback(hObject, eventdata, handles)
% hObject    handle to edt_outputId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_outputId as text
%        str2double(get(hObject,'String')) returns contents of edt_outputId as a double


% --- Executes during object creation, after setting all properties.
function edt_outputId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_outputId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edt_subjectId_Callback(hObject, eventdata, handles)
% hObject    handle to edt_subjectId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_subjectId as text
%        str2double(get(hObject,'String')) returns contents of edt_subjectId as a double


% --- Executes when selected object is changed in bg_listLength.
function bg_listLength_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bg_listLength 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if hObject == handles.rb_list20
    set(handles.pop_listNr, 'String', num2cell(1:45));
elseif hObject == handles.rb_list30
    % must ensure item selected in pop_listNr is valid after change (only 40 lists with 30 sentences)
    set(handles.pop_listNr, 'Value', min(40,get(handles.pop_listNr, 'Value')));
    set(handles.pop_listNr, 'String', num2cell(1:40));
else
    warning('Unhandled selection in button group bg_listLength');
end



function edt_outputFs_Callback(hObject, eventdata, handles)
% hObject    handle to edt_outputFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_outputFs as text
%        str2double(get(hObject,'String')) returns contents of edt_outputFs as a double


% --- Executes during object creation, after setting all properties.
function edt_outputFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_outputFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_device.
function pop_device_Callback(hObject, eventdata, handles)
% hObject    handle to pop_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_device contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_device


% --- Executes during object creation, after setting all properties.
function pop_device_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_playbackFs_Callback(hObject, eventdata, handles)
% hObject    handle to edt_playbackFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_playbackFs as text
%        str2double(get(hObject,'String')) returns contents of edt_playbackFs as a double


% --- Executes during object creation, after setting all properties.
function edt_playbackFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_playbackFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_outputId2_Callback(hObject, eventdata, handles)
% hObject    handle to edt_outputId2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_outputId2 as text
%        str2double(get(hObject,'String')) returns contents of edt_outputId2 as a double


% --- Executes during object creation, after setting all properties.
function edt_outputId2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_outputId2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edt_outputFs2_Callback(hObject, eventdata, handles)
% hObject    handle to edt_outputFs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_outputFs2 as text
%        str2double(get(hObject,'String')) returns contents of edt_outputFs2 as a double


% --- Executes during object creation, after setting all properties.
function edt_outputFs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_outputFs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_monitor.
function cb_monitor_Callback(hObject, eventdata, handles)
% hObject    handle to cb_monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_monitor


% --- Executes on button press in cb_randSeekNoise.
function cb_randSeekNoise_Callback(hObject, eventdata, handles)
% hObject    handle to cb_randSeekNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_randSeekNoise



function edt_olsaTarget_Callback(hObject, eventdata, handles)
% hObject    handle to edt_olsaTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_olsaTarget as text
%        str2double(get(hObject,'String')) returns contents of edt_olsaTarget as a double


% --- Executes during object creation, after setting all properties.
function edt_olsaTarget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_olsaTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_minStep_Callback(hObject, eventdata, handles)
% hObject    handle to edt_minStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_minStep as text
%        str2double(get(hObject,'String')) returns contents of edt_minStep as a double


% --- Executes during object creation, after setting all properties.
function edt_minStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_minStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_nomSpeechLvl.
function cb_nomSpeechLvl_Callback(hObject, eventdata, handles)
% hObject    handle to cb_nomSpeechLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_nomSpeechLvl


% --- Executes on button press in rb_s0nmin90.
function rb_s0nmin90_Callback(hObject, eventdata, handles)
% hObject    handle to rb_s0nmin90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_s0nmin90


% --- Executes on button press in rb_s0nplus70.
function rb_s0nplus70_Callback(hObject, eventdata, handles)
% hObject    handle to rb_s0nplus70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_s0nplus70

% --- Executes on button press in rb_s0nmin70.
function rb_s0nmin70_Callback(hObject, eventdata, handles)
% hObject    handle to rb_s0nmin70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_s0nmin70


% --- Executes on button press in cb_vocoder_L.
function cb_vocoder_L_Callback(hObject, eventdata, handles)
% hObject    handle to cb_vocoder_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_vocoder_L


% --- Executes on button press in cb_vocoder_R.
function cb_vocoder_R_Callback(hObject, eventdata, handles)
% hObject    handle to cb_vocoder_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_vocoder_R


% --- Executes on button press in cb_speech_L.
function cb_speech_L_Callback(hObject, eventdata, handles)
% hObject    handle to cb_speech_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_speech_L


% --- Executes on button press in cb_speech_R.
function cb_speech_R_Callback(hObject, eventdata, handles)
% hObject    handle to cb_speech_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_speech_R


% --- Executes on button press in pupil_dilation.
function pupil_dilation_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_dilation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_dilation

% --- Executes on button press in pupilcheck.
function pupilcheck_Callback(hObject, eventdata, handles)
% hObject    handle to pupilcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupilcheck

% button_state = get(hObject, 'Value');

hFig = matrix_gui_main;
h = guihandles(hFig);

[filename, pathname] = uigetfile({'*.mat'},'Load Data');
if pathname==0
	return
end

load([pathname filename]);
%assignin('base','rgb50',num2str(data.rgb50,0))

set(h.edt_subjectId, 'String',data.subject)
disp('Here is the value for rgb directly from the GUI')
disp(data.rgb50)
set(h.rgb, 'Value', data.rgb50)
%set(h.edt_rgb, 'String', data.rgb50)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% code from sebas script
% --- Executes on button press in openset.
function openset_Callback(hObject, eventdata, handles)
% hObject    handle to openset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of openset



function edt_rgb_Callback(hObject, eventdata, handles)
% hObject    handle to edt_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_rgb as text
%        str2double(get(hObject,'String')) returns contents of edt_rgb as a double


% --- Executes during object creation, after setting all properties.
function edt_rgb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of code from seba


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pupilcheck.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pupilcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edt_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_rgb as text
%        str2double(get(hObject,'String')) returns contents of edt_rgb as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_noise_R.
function cb_noise_R_Callback(hObject, eventdata, handles)
% hObject    handle to cb_noise_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_noise_R


% --- Executes on button press in cb_noise_L.
function cb_noise_L_Callback(hObject, eventdata, handles)
% hObject    handle to cb_noise_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_noise_L


% --- Executes during object creation, after setting all properties.
function fig_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fig_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
