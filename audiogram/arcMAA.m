function varargout = arcMAA(varargin)
% ARCMAA MATLAB code for arcMAA.fig
%      ARCMAA, by itself, creates a new ARCMAA or raises the existing
%      singleton*.
%
%      H = ARCMAA returns the handle to a new ARCMAA or the handle to
%      the existing singleton*.
%
%      ARCMAA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARCMAA.M with the given input arguments.
%
%      ARCMAA('Property','Value',...) creates a new ARCMAA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arcMAA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arcMAA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arcMAA

% Last Modified by GUIDE v2.5 12-Feb-2016 16:13:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arcMAA_OpeningFcn, ...
                   'gui_OutputFcn',  @arcMAA_OutputFcn, ...
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


% --- Executes just before arcMAA is made visible.
function arcMAA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arcMAA (see VARARGIN)


%% Get all important trial and configuration info
handles			= getcfgdefaults(handles);
sethandles(handles);
handles			= gethandles(handles);

%% Initialize TDT
handles			= tdt_init(handles); % Turning on, Loading, and Running Circuits

%% Start experiment
handles.data	= struct([]);
handles			= setupShow(handles);


% Choose default command line output for arcMAA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arcMAA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arcMAA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function sethandles(handles)
set(handles.popupmenu_expinitials,'String',handles.cfg.expinitials);
set(handles.edit_date,'String',handles.cfg.date);

function handles = gethandles(handles)
cfg					= handles.cfg;
expIdx				= get(handles.popupmenu_expinitials,'Value');
expInitials			= get(handles.popupmenu_expinitials,'String');
cfg.expInitials		= expInitials{expIdx};
cfg.subjectid		= get(handles.edit_subjectid,'String');
cfg.datestring		= get(handles.edit_date,'String');
cfg.block			= get(handles.edit_block,'String');
cfg.fname			= [cfg.expInitials '-' sprintf('%04u',str2double(cfg.subjectid)) '-' cfg.datestring '-' sprintf('%04u',str2double(cfg.block)) '.sphere']; % file name
cfg.dname			= [cfg.fpath cfg.expInitials filesep cfg.expInitials '-' sprintf('%04u',str2double(cfg.subjectid)) '-' cfg.datestring]; % directory name
cfg.expdir			= [cfg.fpath cfg.expInitials filesep 'EXP' filesep]; % exp directory name
cfg.snddir			= [cfg.fpath cfg.expInitials filesep 'SND' filesep]; % wav directory name
expIdx				= get(handles.popupmenu_expinitials,'Value');

% get exp and cfg files
str					= [cfg.expdir filesep '*.exp'];
d					= dir(str); % default exp folder
cfg.expfiles		= {d.name};
if isempty(cfg.expfiles)
	
end
% set(handles.popupmenu_exp,'String',cfg.expfiles)
% expfileIdx			= get(handles.popupmenu_exp,'Value');
% cfg.expfname		= cfg.expfiles{expfileIdx};

% d					= dir([cfg.expdir filesep '*.cfg']);
% cfg.cfgfiles		= {d.name};
% set(handles.popupmenu_cfg,'String',cfg.cfgfiles);
% cfgfileIdx			= get(handles.popupmenu_cfg,'Value');
% cfg.cfgfname		= cfg.cfgfiles{cfgfileIdx};
handles.cfg			= cfg;

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

%% standard variables

if ~isfield(cfg,'RZ6_1circuit')
	% RP2_1circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\hrtf_2Channels_rp2.rco';
	% RP2_1circuit = 'C:\Gitlab\Marc\experiment\RPvdsEx\GvB_GWN_Single_Sound.rcx';
	cfg.RZ6_1circuit = which('\MAA_RZ6.rcx');
end

%% short names of TDT circuits
[~,cfg.sRZ6_1circuit]	= fileparts(cfg.RZ6_1circuit);
				
function handles = tdt_init(handles)
% TDTINIT
%
% Initializes TDT system, checks and opens graphical monitor
%
% See also TDT_GLOBALS, TDT_MONITOR


t = tic;

%% Active X Control/Objects
% handles.cfg.RZ6_1circuit = 
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RZ6_1, err(2)]		= RZ6(1,handles.cfg.RZ6_1circuit); % Real-time processor 1
% [RA16_1, err(3)]	= RA16(1,handles.cfg.RA16_1circuit); % Real-time acquisition
for muxIdx = 1:3
	MUX(RZ6_1,muxIdx);
end


%% TDT status
handles.cfg.RZ6_1Status	= RZ6_1.GetStatus;
% handles.cfg.RP2_1Status		= RP2_1.GetStatus;
% handles.cfg.RP2_2Status		= RP2_2.GetStatus;
handles					= tdt_monitor(handles);
% toc(t)
% 
% %% Configuration
% handles.cfg.RP2_1	= RP2_1;
% handles.cfg.RP2_2	= RP2_2;
handles.cfg.RZ6_1	= RZ6_1;
handles.cfg.zBus	= zBus;
% handles.cfg.PA5_1	= PA5_1;
% handles.cfg.PA5_2	= PA5_2;
% handles.cfg.PA5_3	= PA5_3;
% handles.cfg.PA5_4	= PA5_4;

function handles = tdt_monitor(handles)
% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

set(handles.checkbox_RZ6connect,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RZ6load,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RZ6run,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTRun));

function handles = setupShow(handles)
% Modify to your own preference



% Trace: elevation vs azimuth
axes(handles.axes_xy); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
box off
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('azimuth (deg)');
ylabel('elevation (deg)');
title('2D trace');



function runExperiment(handles)

% get trial information
handles = setupShow(handles);
handles		= gethandles(handles);
% handles		= gettrial(handles);
createdir(handles);
% 
% tstart		= tic;
% 
% %%
% trialClean(struct([]),handles.cfg);
handles.cfg.ntrials = 10;
for trlIdx	= 1:handles.cfg.ntrials
	t				= tic;
    
% 	
% 	% Random intertrial interval
% 	r				= rndval(handles.cfg.ITI(1),handles.cfg.ITI(2),1);
% 	handles.trial(trlIdx).ITI		= r; % ms
% 	pause(r/1000);
% 	
% 	%% Trial number
% 	handles.cfg.trial		= trlIdx;
% 	disp(['Trial: ' num2str(trlIdx)])
% 	
% 	%% Trial
% 	stim			= handles.trial(trlIdx).stim;
% 	trialClean(stim,handles.cfg);
% 	
% 	stim			= trialSetup(handles.cfg,stim);
stim = struct([]);
	trialRun(handles.cfg,stim);
    pause(8)
% 	handles.data		= trialSave(handles.cfg,handles.trial,handles.data);
% 	handles				= trialShow(handles);
% 	trialClean(stim,handles.cfg);
% 	
% 	%% End
% 	% save timing
% 	dur		= toc(t);
% 	[~,fname,~]		= fileparts(handles.cfg.fname); % remove extension
% 	fname			= [fname '-' num2str(handles.cfg.trial,'%04u')];
% 	fname			= fcheckext(fname,'sphere');
% 	fname			= fullfile([handles.cfg.dname filesep 'trial' filesep],fname);
% 	save(fname,'dur','-append');
% 
	toc(t)
	disp('==========================================');
end
% 
% 
% handles.cfg.duration	= duration(0,0,toc(tstart)); % add duration of experiment
% handles.cfg				= rmfield(handles.cfg,{'hcurtar','hcurdat'}); % remove figure handles from cfg file
% 
% %% Finish
% % data	= handles.data;
% % cfg		= handles.cfg;
% % trial	= handles.trial;
% % save(fullfile(cfg.dname,cfg.fname),'data','cfg','trial');
% 
% 
% str = {	'DATA ARE SAVED PER TRIAL. SEE ALSO SPHERETRIAL2COMPLETE',...
% 	'______________________________________________________',...
% 	' ',...
% 	['Experiment was completed in: ' char(handles.cfg.duration) ' hours:minutes:seconds'],...
% 	['Files were saved to ' handles.cfg.dname filesep 'trial'],...
% 	};
% str = char(str);
% warndlg(str);
% cd(handles.cfg.dname)
% 
% endBlock(handles.cfg);

function trialRun(cfg,stim)
% RUNTRIAL(ZBUS,RA16)

%% Run zBus
cfg.zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
cfg.zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

% %% Button Press
% disp('Waiting for RA16 button press/sound/led/acquisition');
% % pause(.1); % short break
% while ~cfg.RA16_1.GetTagVal('Wait')
% 	% 	disp('waiting')
% 	% do nothing
% end
% 
% %% Sound play
% % TODO: correct waiting/state machine
% if strcmp(cfg.sRP2_1circuit,'RP2_WAV')
% 	disp('Waiting for RP2 sound');
% 	selsnd = strcmpi({stim.modality},'sound');
% 	if any(selsnd)
% 		snd		= stim(selsnd);
% 		nsnd	= numel(snd);
% 		for sndIdx = 1:nsnd
% 			sndsetup	= cfg.lookup(snd(sndIdx).Z+1,2:4);
% 			RP2out		= num2str(cfg.mux2rp2(sndsetup(2)));
% 			switch sndsetup(1)
% 				case 1
% 					cfg.RP2_1.GetTagVal(['bufferPos' RP2out])
% 					while cfg.RP2_1.GetTagVal(['bufferPos' RP2out])
% 						% do nothing
% 					end
% 				case 2
% 					cfg.RP2_2.GetTagVal(['bufferPos' RP2out])
% 					while cfg.RP2_2.GetTagVal(['bufferPos' RP2out])
% 						% do nothing
% 					end
% 			end
% 		end
% 	end
% 	
% end
% 
% %% Data acquisition
% disp('Waiting for data acquisition');
% % pause(.1); % short break
% while cfg.RA16_1.GetTagVal('Active')
% 	% do nothing
% end


function handles = createdir(handles)
% [FNAME,DNAME] = GETFNAME(FPATH)
%
% Input dialog for data file name
%
% See also FEXISTDLG

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com

%% Initialization

cfg		= handles.cfg;

if ~exist(cfg.dname,'dir') % create directory
	mkdir(cfg.dname);
	mkdir([cfg.dname filesep 'trial']);
	str = ['Creating directory: ' cfg.dname];
	disp(str);
else
	str = {['   Directory ''' cfg.dname ''' already exists.'];'   Data file will be added to this folder.'};
	disp(char(str));
	mkdir([cfg.dname filesep 'trial']);
end

handles.cfg = cfg;

% --- Executes on button press in checkbox_RZ6connect.
function checkbox_RZ6connect_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RZ6connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RZ6connect


% --- Executes on button press in checkbox_RZ6load.
function checkbox_RZ6load_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RZ6load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RZ6load


% --- Executes on button press in checkbox_RZ6run.
function checkbox_RZ6run_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RZ6run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RZ6run


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
runExperiment(handles)


function edit_subjectid_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subjectid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subjectid as text
%        str2double(get(hObject,'String')) returns contents of edit_subjectid as a double


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
