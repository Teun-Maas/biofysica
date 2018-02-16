function varargout = sphereMinor(varargin)
% SPHEREMINOR MATLAB code for sphereMinor.fig
%      SPHEREMINOR, by itself, creates a new SPHEREMINOR or raises the existing
%      singleton*.
%
%      H = SPHEREMINOR returns the handle to a new SPHEREMINOR or the handle to
%      the existing singleton*.
%
%      SPHEREMINOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPHEREMINOR.M with the given input arguments.
%
%      SPHEREMINOR('Property','Value',...) creates a new SPHEREMINOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sphereMinor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sphereMinor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sphereMinor

% Last Modified by GUIDE v2.5 16-Feb-2018 14:45:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @sphereMinor_OpeningFcn, ...
	'gui_OutputFcn',  @sphereMinor_OutputFcn, ...
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

% --- Executes just before sphereMinor is made visible.
function sphereMinor_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sphereMinor (see VARARGIN)


%% Get all important trial and configuration info
handles			= getcfgdefaults(handles);
sethandles(handles);
handles			= gethandles(handles);

%% Initialize TDT
handles			= tdt_init(handles); % Turning on, Loading, and Running Circuits

%% Start experiment
handles.data	= struct([]);
handles			= setupShow(handles);

% Choose default command line output for sphereMinor
handles.output	= hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sphereMinor wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = sphereMinor_OutputFcn(hObject, eventdata, handles)
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
runExperiment(handles)

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


function handles = getcfgdefaults(handles)
%% Default parameters
cfg.fpath			= ['C:' filesep 'DATA' filesep]; % default data folder
cfg.date			= date;

cd(cfg.fpath); % default data folder
d					= dir;
d					= d([d.isdir]);
d					= d(3:end);
cfg.expinitials		= {d.name};

formatOut	= 'yy-mm-dd';
d			= datestr(now,formatOut);
cfg.date	= d;

%% TDT and PLC defaults
cfg = tdt_globals(cfg);		% TDT defaults, in cfg structure, includes circuit names for [RA16_1,RA16_2,RX6,RP2_1,RP2_2]

%% TDT and PLC defaults
cfg.nleds				= 8; % maximum number of LED configurations on PLC
cfg.calfile				= which('sphere.net');
if isempty(cfg.calfile)
	cfg.calfile				= which('defaultsphere.net');
end
%% Convert experimental trial

cfg						= spherelookup(cfg); % lookup structure

%% Data filter coefficients
cfg.Fcutoff             = 80;       % Cutoff Frequency  (Hz)
cfg.Order                   = 50;

cfg.lpFilt = designfilt('lowpassfir', 'FilterOrder', cfg.Order, 'CutoffFrequency', ...
                    cfg.Fcutoff, 'SampleRate', cfg.medusaFs, 'Window', 'hamming');
				
%% Led colours
cfg.ledcolours = {'g','r'}; % green and red

%% handles
handles.cfg		= cfg;



function cfg = tdt_globals(cfg)
% TDT_GLOBALS
%
% "Globals" for TDT system:
% - Circuit filenames for RA16 and RP2s
%
% These circuit filenames are searched on the Matlab path.

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com

%% TDT flags
cfg.statTDTConnect = 1;
cfg.statTDTLoad    = 2;
cfg.statTDTRun     = 3;


%% Fixed RA16.rcx
cfg.RZ6Fs	=  48828.125000; % Hz
cfg.dataidx		= {'Data_1' 'Data_2' 'Data_3' 'Data_4' 'Data_5' 'Data_6' 'Data_7' 'Data_8'}; % names of Data sources

%% RP2 and Muxes
% cfg.mux2rp2				= [1 2 2 1]; % which RP2 channel belongs to which MUX?
cfg.mux2rp2				= [1 1 2 2]; % which RP2 channel belongs to which MUX?

cfg.recdataidx		= {'recData_1' 'recData_2'}; % names of Data sources

%% Fixed setup (SA1 & PA5 & RP2.rcx)
cfg.maxsndlevel		= 75; % for GWN amplitude 1 in RP2.rcx

%% Hardware "defaults" - needs to be checked every time
cfg.SA1gain			= 0;
cfg.remmel.X.gain	= 500;
cfg.remmel.X.offset	= 500;
cfg.remmel.X.invert	= 0;
cfg.remmel.Y.gain	= 500;
cfg.remmel.Y.offset	= 500;
cfg.remmel.Y.invert	= 0;
cfg.remmel.Z.gain	= 500;
cfg.remmel.Z.offset	= 500;
cfg.remmel.Z.invert	= 0;

%% standard variables
if ~isfield(cfg,'RA16_1circuit')
	% RA16_1circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\Remmel_8Channels_ra16_1.rco';
	% RA16_1circuit = 'C:\MATLAB\experiment\RPvdsEx\GvB_Remmel_8Channels_light_version.rcx';
	cfg.RA16_1circuit = which('sphere_RA16.rcx');
end
if ~isfield(cfg,'RP2_1circuit')
	% RP2_1circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\hrtf_2Channels_rp2.rco';
	% RP2_1circuit = 'C:\Gitlab\Marc\experiment\RPvdsEx\GvB_GWN_Single_Sound.rcx';
	cfg.RP2_1circuit = which('sphere_RP2_WAV.rcx');
end
if ~isfield(cfg,'RP2_2circuit')
	% RP2_2circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\hrtf_2Channels_rp2.rco';
	% RP2_2circuit = 'C:\Gitlab\Marc\experiment\RPvdsEx\GvB_GWN_Single_Sound.rcx';
	cfg.RP2_2circuit = which('sphere_RP2_WAV.rcx');
end

%% short names of TDT circuits
[~,cfg.sRA16_1circuit]	= fileparts(cfg.RA16_1circuit);
[~,cfg.sRP2_1circuit]	= fileparts(cfg.RP2_1circuit);
[~,cfg.sRP2_2circuit]	= fileparts(cfg.RP2_2circuit);


function cfg = spherelookup(cfg)
% CFG = SPHERELOOKUP
%
% Get lookup table for LEDs and speakers in the Sphere setup
%
% CFG contains field
%	- lookup: % [PLC-Chan# RP2# MUX# BIT# AZ EL]

if nargin<1
	cfg = [];
end
%% LED connections
% [Azimuth (deg) Elevation (deg) LED #]
% azimuths and elevations were measured by Sebastian Ausili in June 2015
cfg					= sphere_speakerpositions(cfg);
cfg.nstimchan		= 2^7; % number of PLC and MUX channels
cfg.nspeakers		= 112; % actual number of speakers
cfg.stimchan		= (1:cfg.nstimchan)-1;

%% Add missing channels
n					= size(cfg.lookup,1);
cfg.lookup			= [cfg.lookup; NaN(cfg.nstimchan-n,3)]; % add missing channel data as NaNs
sel					= ismember(cfg.stimchan,cfg.lookup(:,3)); % lookup existing channels
cfg.missingchan		= cfg.stimchan(~sel); % get missing channels
sel					= isnan(cfg.lookup(:,3)); % search for missing channels in lookup-table
cfg.lookup(sel,3)	= cfg.missingchan; % put missing channel numbers in lookup-table
cfg.lookup			= sortrows(cfg.lookup,[3 1 2]); % sort lookup-table by channel number

%% Speakers
% [Channel# RP2# MUX# BIT#]
%
% RP2ind   = [1 2];
% MUXind   = 1:4;
% MUXbit   = 1:16;
% ChanNo  = 0:127;
% LookUp = NaN(128.4);
% LookUp(:.1) = ChanNo;
% LookUp(1:64.2) = RP2ind(1);
% LookUp(65:128.2) = RP2ind(2);
% LookUp(1:16.3) = MUXind(1);
% LookUp(17:32.3) = MUXind(2);
% LookUp(33:48.3) = MUXind(3);
% LookUp(49:64.3) = MUXind(4);
% LookUp(65:80.3) = MUXind(1);
% LookUp(81:96.3) = MUXind(2);
% LookUp(97:112.3) = MUXind(3);
% LookUp(113:128.3) = MUXind(4);
% LookUp(:.4) = [MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit]';
cfg.nRP2			= 2; % 2 RP2s
cfg.nMUX			= 4; % 4 multiplexers
cfg.nMUXbit			= 16; % 16 bits per MUX
cfg.RP2ind			= 1:cfg.nRP2;
cfg.MUXind			= 1:cfg.nMUX;
cfg.MUXbit			= 1:cfg.nMUXbit;
LookUp(:,1)			= cfg.stimchan; %
for ii				= 1:cfg.nRP2
	idx				= (1:cfg.nMUX*cfg.nMUXbit)+(ii-1)*cfg.nMUX*cfg.nMUXbit;
	LookUp(idx,2)	= cfg.RP2ind(ii);
	for jj			= 1:cfg.nMUX
		idx				= (1:cfg.nMUXbit)+(jj-1)*cfg.nMUXbit+(ii-1)*cfg.nMUX*cfg.nMUXbit;
		LookUp(idx,3)	= cfg.MUXind(jj);
	end
end
LookUp(:,4)			= repmat(cfg.MUXbit,1,cfg.nRP2*cfg.nMUX);

%% Combine
% [PLC-Chan# RP2# MUX# BIT# AZ EL]
cfg.lookup		= [cfg.lookup(:,3) LookUp(:,2:4) cfg.lookup(:,1)  cfg.lookup(:,2)];
cfg.lookuplabel = {'Channel' 'RP2' 'MUX' 'Bit' 'Azimuth' 'Elevation'};

%% Channel 79 is missing
% quick fix, replace with Channel 31
% Date: 20-8-2015
cfg.lookup(32,5:6) = cfg.lookup(80,5:6);
cfg.lookup(80,5:6) = NaN;

%% Extra infrared
cfg.lookup(128,5:6) = [700,700];
% keyboard
%% Tannoy speakers
% Januray 8
TannoyLocations				= [-70:10:-10 10:10:70];
TannoyIDs					= [26:30 58:63  124:126];
cfg.TannoyIDs				= TannoyIDs;

cfg.lookup(TannoyIDs+1,5)	= TannoyLocations; % azimuth
cfg.lookup(TannoyIDs+1,6)	= -7.5; % elevation


%% Interpolant
x		= cfg.lookup(:,5);
y		= cfg.lookup(:,6);
z		= cfg.lookup(:,1);
sel		= ~isnan(x);
x		= x(sel);
y		= y(sel);
z		= z(sel);
cfg.interpolant		= scatteredInterpolant(x,y,z,'nearest');
