function varargout = psychoILD(varargin)
% PSYCHOILD MATLAB code for psychoILD.fig
%      PSYCHOILD, by itself, creates a new PSYCHOILD or raises the existing
%      singleton*.
%
%      H = PSYCHOILD returns the handle to a new PSYCHOILD or the handle to
%      the existing singleton*.
%
%      PSYCHOILD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSYCHOILD.M with the given input arguments.
%
%      PSYCHOILD('Property','Value',...) creates a new PSYCHOILD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before psychoILD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to psychoILD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help psychoILD

% Last Modified by GUIDE v2.5 04-Apr-2017 10:26:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @psychoILD_OpeningFcn, ...
	'gui_OutputFcn',  @psychoILD_OutputFcn, ...
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


% --- Executes just before psychoILD is made visible.
function psychoILD_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to psychoILD (see VARARGIN)


%% Get all important trial and configuration info
handles			= getcfgdefaults(handles);
sethandles(handles);
handles			= gethandles(handles);

%% Initialize TDT
handles			= tdt_init(handles); % Turning on, Loading, and Running Circuits

%% Start experiment
handles.data	= struct([]);
handles			= setupShow(handles);

% Choose default command line output for psychoILD
handles.output	= hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes psychoILD wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = psychoILD_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function runExperiment(handles)

% get trial information
handles     = setupShow(handles);
handles		= gethandles(handles);
handles		= gettrial(handles);
createdir(handles);

tstart		= tic;

%%
trialClean(struct([]),handles.cfg);
for trlIdx	= 1:handles.cfg.ntrials
	t				= tic;
	
	% Random intertrial interval
% 	r				= rndval(handles.cfg.ITI(1),handles.cfg.ITI(2),1);
% 	handles.trial(trlIdx).ITI		= r; % ms
% 	pause(r/1000);
	
	%% Trial number
	handles.cfg.trial		= trlIdx;
	disp(['Trial: ' num2str(trlIdx)])
	
	%% Trial
	stim.ild             = handles.trial.ild(trlIdx);
	stim.freq_left		= handles.trial.freq_left(trlIdx);
	stim.freq_right		= handles.trial.freq_right(trlIdx);
	trialClean(stim,handles.cfg);
	
    trialSetup(handles.cfg,stim);
    trialRun(handles.cfg,stim);
	trialSave(handles.cfg,handles.trial);
% 	handles				= trialShow(handles);
	trialClean(struct([]),handles.cfg);
	
	%% End
	% save timing
	dur		= toc(t);
	[~,fname,~]		= fileparts(handles.cfg.fname); % remove extension
	fname			= [fname '-' num2str(handles.cfg.trial,'%04u')]; %#ok<*AGROW>
	fname			= fcheckext(fname,'sphere');
	fname			= fullfile([handles.cfg.dname filesep 'trial' filesep],fname);
% 	save(fname,'dur','-append');

	toc(t)
	disp('==========================================');
end


handles.cfg.duration	= duration(0,0,toc(tstart)); % add duration of experiment
% handles.cfg				= rmfield(handles.cfg,{'hcurtar','hcurdat'}); % remove figure handles from cfg file

%% Finish
% data	= handles.data;
% cfg		= handles.cfg;
% trial	= handles.trial;
% save(fullfile(cfg.dname,cfg.fname),'data','cfg','trial');


str = {	'DATA ARE SAVED PER TRIAL. SEE ALSO SPHERETRIAL2COMPLETE',...
	'______________________________________________________',...
	' ',...
	['Experiment was completed in: ' char(handles.cfg.duration) ' hours:minutes:seconds'],...
	['Files were saved to ' handles.cfg.dname filesep 'trial'],...
	};
str = char(str);
% warndlg(str);
disp(str);
cd(handles.cfg.dname)

% endBlock(handles.cfg);

function trialSetup(cfg,stim)

setSound(stim,cfg)

function setSound(snd,cfg)
% SETSOUND(SND,CFG,RP2STR)
%
% Set all parameters for sound presentation

cfg.RZ6_1.SetTagVal('FC1',snd.freq_left);
cfg.RZ6_1.SetTagVal('FC2',snd.freq_right);

attA = snd.ild/2+50;
attB = -snd.ild/2+50;

cfg.RZ6_1.SetTagVal('attA',attA);
cfg.RZ6_1.SetTagVal('attB',attB);

function endBlock(cfg)


snd.X			= 0;
snd.Y			= 0;
snd.channel		= [];
snd.detect		= [];
snd.event		= [];
snd.intensity	= 60;
snd.modality	= 'sound';
snd.offdelay	= 0;
snd.offevent	= 0;
snd.ondelay		= 0;
snd.onevent		= 0;
snd.matfile		= 'rehandel.mat';
snd.azimuth		= 3.5084e-15;
snd.elevation	= -3.5084e-15;
snd.Z			= 10;
snd.ledhandle	= [];

% load handel;
% y = resample(y,48828,Fs);
% snd = 
% fname		= fullfile(cfg.snddir,snd.matfile);
% if ~exist(fname,'file')
% 	disp('You do not have Handel installed.');
% 	disp('1000 monkeys are doing your work for you, and are now writing Handel.');
% 	load handel;
% 	snd = resample(y,48828,Fs);
% 	save(fname,'snd');
% end
% % setSound(snd,cfg,'RP2_1');
% stim = trialSetup(cfg,snd);
% 
% trialRun(cfg,snd);

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
set(handles.popupmenu_exp,'String',cfg.expfiles)
expfileIdx			= get(handles.popupmenu_exp,'Value');

cfg.expfname		= cfg.expfiles{expfileIdx};

d					= dir([cfg.expdir filesep '*.cfg']);
cfg.cfgfiles		= {d.name};
set(handles.popupmenu_cfg,'String',cfg.cfgfiles);
cfgfileIdx			= get(handles.popupmenu_cfg,'Value');
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

function handles = gettrial(handles)
cfg			= handles.cfg;

% %% Data file name dialog and file check
% handles					= getfname(handles); % datafile name dialog
cfg						= fexistdlg(handles.cfg); % check datafile name
% cfg.dname				= [cfg.dname filesep];
% if isempty(cfg.fname);
% 	return % quit if there is no file name
% end
%
% %% Exp-file
% % Dialog to choose the exp-file and load experimental parameters
[trial,cfg]				= readexp(cfg); % load experimental trial parameters and configuration

%% handles
handles.trial	= trial;
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


%% Fixed RZ6.rcx sample rate
cfg.Fs	= 48828.125; % Hz


%% Fixed setup (SA1 & PA5 & RP2.rcx)
cfg.maxsndlevel		= 75; % for GWN amplitude 1 in RP2.rcx


%% standard variables
if ~isfield(cfg,'RZ6circuit')
	% RA16_1circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\Remmel_8Channels_ra16_1.rco';
	% RA16_1circuit = 'C:\MATLAB\experiment\RPvdsEx\GvB_Remmel_8Channels_light_version.rcx';
	cfg.RZ6_1circuit = which('psychoILD.rcx');
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
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RZ6_1, err(2)]	= RZ6(1,handles.cfg.RZ6_1circuit); % Real-time acquisition

%% TDT status
handles.cfg.RZ6_1Status	= RZ6_1.GetStatus;
handles					= tdt_monitor(handles);
toc(t)

%% Configuration
handles.cfg.RZ6_1	= RZ6_1;
handles.cfg.zBus	= zBus;

function handles = tdt_monitor(handles)
% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

set(handles.checkbox_RZ6connect,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RZ6load,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RZ6run,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTRun));
set(handles.text_RZ6,'String',['RZ6 #1: ' handles.cfg.sRZ6_1circuit '.rcx']);

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

function cfg = fexistdlg(cfg)
% [FNAME,DNAME] = FEXISTDLG(FNAME,DNAME)
%
% Check whether the current data file name does not already exist. If so,
% prompts experimenter to continue, stop, or rename.
%
% See also GETFNAME

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com
fileexist		= exist([cfg.dname filesep cfg.fname],'file');
while fileexist
	choice1 = 'Continue/overwrite';
	choice2 = 'Stop';
	choice		= questdlg(['File ' cfg.fname ' already exists. What dou you want to do?'],'Check filename',choice1,choice2,choice2);
	switch choice
		case(choice1)
			disp('Continuing...');
			fileexist = false;
		case(choice2)
			% 			disp('Quitting...');
			cfg.fname = [];
			cfg.dname = [];
			return
	end
end

function handles = setupShow(handles)
% Modify to your own preference



% Trace: elevation vs azimuth
axes(handles.axes_psi); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
box off
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('azimuth (deg)');
ylabel('elevation (deg)');
title('2D trace');


axes(handles.axes_button); %#ok<*NASGU>
cla
hold on
box off
set(gca,'YTick',-90:30:90,'TickDir','out');
xlabel('Time (s)');
ylabel('Position (deg)');
ylim([-100 +100]);

function handles = trialShow(handles)
% Modify to your own preference

% need for automatic calibration
% need for automatic detection


%% calibrated traces
% history in gray
% current black bold
handles.stim	= handles.trial(handles.cfg.trial).stim;
handles			= stimShow(handles);
handles			= datShow(handles);

%% Linear Regression
% subplot(132)
% subplot(133)

%% Start drawing
drawnow;

function handles = stimShow(handles)

axes(handles.axes_psi); %#ok<*NASGU>
if isfield(handles.cfg,'hcurtar')
	n = numel(handles.cfg.hcurtar);
	for ii = 1:n
		if isfield(get(handles.cfg.hcurtar(ii)),'MarkerSize')
			set(handles.cfg.hcurtar(ii),'MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);
		end
	end
end
for stmIdx = 1:numel(handles.stim)
	switch upper(handles.stim(stmIdx).modality)
		case {'LED','SKY'}
			handles.cfg.hcurtar(stmIdx) = plot(handles.stim(stmIdx).azimuth,handles.stim(stmIdx).elevation,'r*','MarkerSize',8);
		case {'SOUND','SND','SND1','SND2','SND3','SND4'}
			handles.cfg.hcurtar(stmIdx) = plot(handles.stim(stmIdx).azimuth,handles.stim(stmIdx).elevation,'bo','MarkerSize',8);
	end
end

function handles = datShow(handles)


%% load data
Hchan		= 3; Vchan		= 1; Fchan		= 2; % set channels in correct order
S			= load(handles.cfg.calfile,'-mat');
dat			=  handles.data(1).raw;
H			= dat(:,Hchan); H = H(:);
V			= dat(:,Vchan); V = V(:);
F			= dat(:,Fchan); F = F(:);
DAT			= [H V F]';

%% calibrate
[az,el]		= pa_calib(DAT,S); % calibrate: from Volts 2 Deg
az          = filtfilt(handles.cfg.lpFilt,az); % filter
el          = filtfilt(handles.cfg.lpFilt,el); % filter

t			= (1:numel(az))/handles.cfg.medusaFs; % time (s)
[vel,smv]	= getvel(az',el',handles.cfg.medusaFs,0.01); % get velocities

%% plot
if isfield(handles.cfg,'hcurdat')
	n = numel(handles.cfg.hcurdat);
	for figIdx = 1:n
		if ismember(figIdx,[1:3 5])
			set(handles.cfg.hcurdat(figIdx),'LineWidth',1,'Color',[.7 .7 .7]);
		else
			delete(handles.cfg.hcurdat(figIdx));
		end
	end
end
axes(handles.axes_psi);
handles.cfg.hcurdat(1) = plot(az,el,'k-','lineWidth',2);

axes(handles.axes_button);
handles.cfg.hcurdat(2) = plot(t,az,'b-','lineWidth',2);
hold on
handles.cfg.hcurdat(3) = plot(t,el,'r-','lineWidth',2);
title(handles.cfg.trial);

axes(handles.axes_velocity);
handles.cfg.hcurdat(4) = plot(t,vel,'g-','lineWidth',2);
hold on
handles.cfg.hcurdat(5) = plot(t,smv,'k-','lineWidth',2);

function trialRun(cfg,stim)
% RUNTRIAL(ZBUS,RA16)

%% Run zBus
cfg.zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
cfg.zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

%% Button Press
disp('Waiting for RZ6 button press/sound/led/acquisition');

pause(3)
% pause(.1); % short break
% while ~cfg.RZ6_1.GetTagVal('Wait')
% 	% 	disp('waiting')
% 	% do nothing
% end



% %% Data acquisition
% disp('Waiting for data acquisition');
% % pause(.1); % short break
% while cfg.RZ6_1.GetTagVal('Active')
% 	% do nothing
% end

function button = trialSave(cfg,trial)
% DATA = SAVETRIAL(RA16_1,CFG,LED,DATA)
%
% Save data


%% Read RA16 Event buffer
button				= cfg.RZ6_1.GetTagVal('ButtonPress');
reactiontime		= cfg.RZ6_1.GetTagVal('ReactionTime'); % samples
SF					= cfg.RZ6_1.GetSFreq; % Hz
RT					= reactiontime/SF; % seconds


% %% Save trial data
% % Because an experiment might be stopped at any time, we need to record
% % data trial by trial in order to prevent data loss.
[~,fname,~]		= fileparts(cfg.fname); % remove extension
fname			= [fname '-' num2str(cfg.trial,'%04u')];
% % sndname = fname;
fname			= fcheckext(fname,'ild');
fname			= fullfile([cfg.dname filesep 'trial' filesep],fname);
% 
% %% remove graphics
% if isfield(cfg,'hcurtar') 
% 	cfg				= rmfield(cfg,'hcurtar');
% elseif isfield(cfg,'hcurdat') 
% 	cfg				= rmfield(cfg,'hcurdat');
% else
% 	warning('No graphic handles present to be removed');
% end
% 
stim.ild			= trial.ild(cfg.trial);
stim.freq_left		= trial.freq_left(cfg.trial);
stim.freq_right	= trial.freq_right(cfg.trial); %#ok<*STRNU>

save(fname,'stim','cfg','button','RT');


function trialClean(stim,cfg)
% %% Remove ledhandle if it exists
% nstim = numel(stim);
% for stmIdx = 1:nstim
% 	if isfield(stim(stmIdx),'ledhandle')
% 		if ~isempty(stim(stmIdx).ledhandle)
% 			stim(stmIdx).ledhandle.delete; 	% delete(leds)/switch off light;
% 		end
% 	end
% end
% 
% %% Turn sounds off
% % by switching off bits on inactive PM2Relay muliplexers
% for muxIdx = 1:4
% 	MUX(cfg.RP2_1,muxIdx,0)
% 	MUX(cfg.RP2_2,muxIdx,0)
% end
% 
% % and by setting buffersize to 0
% for rpIdx = 1:2
% 	RPstr = ['RP2_' num2str(rpIdx)];
% 	cfg.RA16_1.SetTagVal(['rp2Enable' num2str(rpIdx)],0);
% 	for chanIdx = 1:2
% 		str		= ['bufferSize' num2str(chanIdx)];
% 		cfg.(RPstr).SetTagVal(str,0);
% 		cfg.(RPstr).SetTagVal(['chanEnable'  num2str(chanIdx)],0);
% 	end
% end

%% Reset of events to unlikely high number
% unlikely = 100;
% for RPidx = 1:2
% 	str = ['eventRP2' num2str(RPidx)];
% 	RA16_1.SetTagVal(str,unlikely);
% end
%
% for LEDidx = 1:8
% 	str = ['eventLED' num2str(LEDidx)];
% 	RA16_1.SetTagVal(str,unlikely);
% end
% RA16_1.SetTagVal('eventWait',unlikely);
% RA16_1.SetTagVal('eventSample',unlikely);

%% Reset
cfg.zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

function [trial,cfg] = readexp(cfg)
% [TRIAL,CFG] = READEXP(FNAME)
%
% Load TRIAL and ConFiGuration parameters from exp-file FNAME
%
% See also GETEXP

% 2015 Marc van Wanrooij

%% Initialization
if nargin<1
	fname = fcheckexist([]);
	cfg = [];
else
	if isstruct(cfg)
		fname = fullfile(cfg.expdir,cfg.expfname);
	else
		fname	= cfg;
		cfg		= [];
	end
end

%% Open
load(fname,'-mat');

trial.ild = ild;
trial.freq_left = freq_left;
trial.freq_right = freq_right;

cfg.ntrials = numel(ild);


function [veltrace,smvtrace,acctrace,smatrace]          = getvel(htrace,vtrace,Fsample,sd)
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% [VEL, SVEL] = GETVEL(HOR,VER,SMOOTHFACTOR,FSAMPLE)
%
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% See also GSMOOTH, PA_SACDET
%
% MarcW 2007

Rx                                      = htrace;
Ry                                      = vtrace;
R                                       = NaN*Rx;
veltrace                                = R;
acctrace                                = R;
smvtrace                                = R;
smatrace                                = R;
for i                                   = 1:size(htrace,2)
	Rx(:,i)                             = gradient(Rx(:,i),1);
	Ry(:,i)                             = gradient(Ry(:,i),1);
	R(:,i)                              = hypot(Rx(:,i),Ry(:,i));
	R(:,i)                              = cumsum(R(:,i));
	
	veltrace(:,i)                       = gradient(R(:,i),1./Fsample);
	smvtrace(:,i)                       = pa_gsmooth(veltrace(:,i),Fsample,sd);
	
	acctrace(:,i)                       = gradient(smvtrace(:,i),1./Fsample);
	smatrace(:,i)                       = pa_gsmooth(acctrace(:,i),Fsample,sd);
end

% --- Executes on button press in checkbox_RZ6connect.
function checkbox_RZ6connect_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
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
