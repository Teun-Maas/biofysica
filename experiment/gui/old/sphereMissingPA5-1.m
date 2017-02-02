function varargout = sphereMissingPA5(varargin)
% SPHEREMISSINGPA5 MATLAB code for sphereMissingPA5.fig
%      SPHEREMISSINGPA5, by itself, creates a new SPHEREMISSINGPA5 or raises the existing
%      singleton*.
%
%      H = SPHEREMISSINGPA5 returns the handle to a new SPHEREMISSINGPA5 or the handle to
%      the existing singleton*.
%
%      SPHEREMISSINGPA5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPHEREMISSINGPA5.M with the given input arguments.
%
%      SPHEREMISSINGPA5('Property','Value',...) creates a new SPHEREMISSINGPA5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sphereMissingPA5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sphereMissingPA5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sphereMissingPA5

% Last Modified by GUIDE v2.5 04-Jan-2016 16:09:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @sphereMissingPA5_OpeningFcn, ...
	'gui_OutputFcn',  @sphereMissingPA5_OutputFcn, ...
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


% --- Executes just before sphereMissingPA5 is made visible.
function sphereMissingPA5_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sphereMissingPA5 (see VARARGIN)


%% Get all important trial and configuration info
handles			= getcfgdefaults(handles);
sethandles(handles);
handles			= gethandles(handles);

%% Initialize TDT
handles			= tdt_init(handles); % Turning on, Loading, and Running Circuits

%% Start experiment
handles.data	= struct([]);
handles			= setupShow(handles);

% Choose default command line output for sphereMissingPA5
handles.output	= hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sphereMissingPA5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = sphereMissingPA5_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function runExperiment(handles)

% get trial information
handles = setupShow(handles);
handles		= gethandles(handles);
handles		= gettrial(handles);
createdir(handles);

tstart		= tic;

%%
trialClean(struct([]),handles.cfg);
for trlIdx	= 1:handles.cfg.ntrials
	t				= tic;
	
	% Random intertrial interval
	r				= rndval(handles.cfg.ITI(1),handles.cfg.ITI(2),1);
	handles.trial(trlIdx).ITI		= r; % ms
	pause(r/1000);
	
	%% Trial number
	handles.cfg.trial		= trlIdx;
	disp(['Trial: ' num2str(trlIdx)])
	
	%% Trial
	stim			= handles.trial(trlIdx).stim;
	trialClean(stim,handles.cfg);
	
	stim			= trialSetup(handles.cfg,stim);
	trialRun(handles.cfg,stim);
	handles.data		= trialSave(handles.cfg,handles.trial,handles.data);
	handles				= trialShow(handles);
	trialClean(stim,handles.cfg);
	
	%% End
	% save timing
	dur		= toc(t);
	[~,fname,~]		= fileparts(handles.cfg.fname); % remove extension
	fname			= [fname '-' num2str(handles.cfg.trial,'%04u')];
	fname			= fcheckext(fname,'sphere');
	fname			= fullfile([handles.cfg.dname filesep 'trial' filesep],fname);
	save(fname,'dur','-append');

	toc(t)
	disp('==========================================');
end


handles.cfg.duration	= duration(0,0,toc(tstart)); % add duration of experiment
handles.cfg				= rmfield(handles.cfg,{'hcurtar','hcurdat'}); % remove figure handles from cfg file

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
warndlg(str);
cd(handles.cfg.dname)

endBlock(handles.cfg);

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
fname		= fullfile(cfg.snddir,snd.matfile);
if ~exist(fname,'file')
	disp('You do not have Handel installed.');
	disp('1000 monkeys are doing your work for you, and are now writing Handel.');
	load handel;
	snd = resample(y,48828,Fs);
	save(fname,'snd');
end
% setSound(snd,cfg,'RP2_1');
stim = trialSetup(cfg,snd);

trialRun(cfg,snd);

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
cfg.cfgfname		= cfg.cfgfiles{cfgfileIdx};
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
				
%% Note
cfg.note	= {'This is a script for a missing PA54',date};

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
%
% %% CFG file
cfg						= readcfg(cfg); % read cfg cfile
cfg.acqdur				= cfg.humanv1.ADC(1).samples / cfg.humanv1.ADC(1).rate * 1000; % TODO: HumanV1/duration of data acquisition (ms)
cfg.nsamples			= round(cfg.acqdur/1000*cfg.medusaFs); % length data acquisition (samples)
cfg.nchan				= 8;


% TODO: code to verify experiment
trial					= sphereZ(trial,cfg);

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


%% Fixed RA16.rcx
cfg.medusaFs	= 6103.515625; % Hz
cfg.dataidx		= {'Data_1' 'Data_2' 'Data_3' 'Data_4' 'Data_5' 'Data_6' 'Data_7' 'Data_8'}; % names of Data sources

%% RP2 and Muxes
cfg.mux2rp2				= [1 2 2 1]; % which RP2 channel belongs to which MUX?
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
[RP2_1, err(2)]		= RP2(1,handles.cfg.RP2_1circuit); % Real-time processor 1
[RP2_2, err(3)]		= RP2(2,handles.cfg.RP2_2circuit); % Real-time processor 2
[RA16_1, err(4)]	= RA16(1,handles.cfg.RA16_1circuit); % Real-time acquisition
[PA5_1, err(5)]		= PA5(1); % Programmable attenuator 1
[PA5_2, err(6)]		= PA5(2); %  Programmable attenuator 2
[PA5_3, err(7)]		= PA5(3); %  Programmable attenuator 3
% [PA5_4, err(8)]		= PA5(4); %  Programmable attenuator 4
for muxIdx = 1:4
	MUX(RP2_1,muxIdx);
	MUX(RP2_2,muxIdx);
end

%% TDT status
handles.cfg.RA16_1Status	= RA16_1.GetStatus;
handles.cfg.RP2_1Status		= RP2_1.GetStatus;
handles.cfg.RP2_2Status		= RP2_2.GetStatus;
handles					= tdt_monitor(handles);
toc(t)

%% Configuration
handles.cfg.RP2_1	= RP2_1;
handles.cfg.RP2_2	= RP2_2;
handles.cfg.RA16_1	= RA16_1;
handles.cfg.zBus	= zBus;
handles.cfg.PA5_1	= PA5_1;
handles.cfg.PA5_2	= PA5_2;
handles.cfg.PA5_3	= PA5_3;
% handles.cfg.PA5_4	= PA5_4;

function handles = tdt_monitor(handles)
% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

set(handles.checkbox_RA16connect,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RA16load,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RA16run,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTRun));

set(handles.checkbox_RP21connect,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RP21load,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RP21run,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTRun));

set(handles.checkbox_RP22connect,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RP22load,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RP22run,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTRun));
set(handles.text_RA16,'String',['RA16 #1: ' handles.cfg.sRA16_1circuit '.rcx']);
set(handles.text_RP21,'String',['RP2 #1: ' handles.cfg.sRP2_1circuit '.rcx']);
set(handles.text_RP22,'String',['RP2 #2: ' handles.cfg.sRP2_2circuit '.rcx']);

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

function trial = sphereZ(trial,cfg)
% TRIAL = SPHEREZ(TRIAL,CFG)
% Convert azimuth, elevation into Z values for sphere
%
% Note that X,Y values do not exist for the SPHERE, only the HOOP
%
% Should be made more accesible
% e.g. not based on structures
%

for trlIdx	= 1:cfg.ntrials
	for stmIdx	= 1:trial(trlIdx).nstim
		X			= trial(trlIdx).stim(stmIdx).X;
		if ~isempty(X)
			ZI										= cfg.interpolant(trial(trlIdx).stim(stmIdx).azimuth,trial(trlIdx).stim(stmIdx).elevation);
			trial(trlIdx).stim(stmIdx).Z			= ZI;
			trial(trlIdx).stim(stmIdx).azimuth		= cfg.lookup(ZI+1,5);
			trial(trlIdx).stim(stmIdx).elevation	= cfg.lookup(ZI+1,6);
			
		end
	end
end

function stim = trialSetup(cfg,stim)
% HLED = RUNTRIAL(RA16)
%
% Set up experimental parameters

%% Set TDT parameters
selled		= strcmpi({stim.modality},'LED') |  strcmpi({stim.modality},'SKY') |  strcmpi({stim.modality},'LASER');
selacq		= strcmpi({stim.modality},'data acquisition');
% seltrg = strcmpi({stim.modality},'trigger');
selsnd		= strcmpi({stim.modality},'sound');
selsndacq	= strcmpi({stim.modality},'sound acquisition');


%% LED
if any(selled)
	led		= stim(selled);
	nled	= numel(led);
	% 	nled = 2
	n		= nled*2; % LEDs need to be turned on and off
	s		= ledpattern(n);
	cnt		= 0;
	for ledIdx = 1:nled
		% TDT RA16
		% Set timing information on LEDs
		% Note that in RA16 circuit, event 1 = start of experiment
		str1 = ['eventLED' num2str(2*ledIdx-1)];
		str2 = ['eventLED' num2str(2*ledIdx)];
		cfg.RA16_1.SetTagVal(str1,led(ledIdx).onevent+1);
		cfg.RA16_1.SetTagVal(str2,led(ledIdx).offevent+1);
		str1 = ['delayLED' num2str(2*ledIdx-1)];
		str2 = ['delayLED' num2str(2*ledIdx)];
		cfg.RA16_1.SetTagVal(str1,led(ledIdx).ondelay+1);
		cfg.RA16_1.SetTagVal(str2,led(ledIdx).offdelay+1);
		
		% PLC
		for ii	= 1:2:n
			cnt = cnt+1;
			s(cnt).set(led(ledIdx).Z,'g');
			s(cnt).intensity('g',round(led(ledIdx).intensity/6)); % hoop: range 0-255, sphere range 1-50
		end
	end
	stim(find(selled,1)).ledhandle = ledcontroller;
	stim(find(selled,1)).ledhandle.write(s);
end

%% Acquisition
if any(selacq)
	acq	= stim(selacq);
	cfg.RA16_1.SetTagVal('eventAcq',acq.onevent+1);
	cfg.RA16_1.SetTagVal('delayAcq',acq.ondelay);
	cfg.RA16_1.SetTagVal('acqSamples',cfg.nsamples); % amount of DA samples
end


%% Sound
% 			[RP2tag1,RP2tag2,~,MUXind,MUXbit1,SpeakerChanNo] = GvB_SoundSpeakerLookUp(azrnd(ii),elrnd(ii),RP2_1,RP2_2,LedLookUpTable);
% 			GvB_MUXSet(RP2tag1,RP2tag2,MUXind,MUXbit1,'set');
if any(selsnd)
	snd		= stim(selsnd);
	nsnd	= numel(snd);
	for sndIdx = 1:nsnd
		sndsetup	= cfg.lookup(snd(sndIdx).Z+1,2:4);
		switch sndsetup(1)
			case 1
				maxSamples = setSound(snd(sndIdx),cfg,'RP2_1');
			case 2
				maxSamples = setSound(snd(sndIdx),cfg,'RP2_2');
		end
	end
end

if ~exist('maxSamples','var')
	maxSamples = 0;
end

%% Acquisition
if any(selsndacq)
	sndacq	= stim(selsndacq);
	cfg.RP2_1.SetTagVal('recBufSize',maxSamples+1000); % amount of DA samples
end

%% Wait for?
% This needs some tweaking
% search for latest event with longest offset
% which should also include sampling period and sound although this does not have an
% offevent

e				= [stim.offevent];
d				= [stim.offdelay];
[mxevent,idx]	= max(e);

mxdelay			= max([d(idx) ceil(1000*cfg.nsamples./cfg.medusaFs) ceil(1000*maxSamples/48828.125)]);

%%
cfg.RA16_1.SetTagVal('eventWait',mxevent+1);
cfg.RA16_1.SetTagVal('delayWait',mxdelay);

function maxSamples = setSound(snd,cfg,RP2str)
% SETSOUND(SND,CFG,RP2STR)
%
% Set all parameters for sound presentation

Z			= snd.Z;
atten		= max(cfg.maxsndlevel-snd.intensity,0);

dur			= snd.offdelay-snd.ondelay;
sndsetup	= cfg.lookup(Z+1,2:4);


%% WAV
if strcmp(cfg.sRP2_1circuit,'sphere_RP2_WAV')
	disp('Sounds are obtained from MAT file');
	fname		= fullfile(cfg.snddir,snd.matfile);
	s			= load(fname); % Read sound from mat file
	sig			= transpose(s.snd);
	maxSamples	= length(sig);
	maxSamples	= min([maxSamples,500000]); % maximum set in RP2 circuit
	dur			= round(1000*maxSamples/cfg.medusaFs);
	RP2out		= num2str(cfg.mux2rp2(sndsetup(2))); % from MUX to RP2 channel
	if strcmp(RP2str,'RP2_1')
		RP2out = '2';
	end
	cfg.(RP2str).SetTagVal(['bufferSize' RP2out],maxSamples);
	cfg.(RP2str).WriteTagV(['wavData' RP2out], 0,sig(1:maxSamples));
	cfg.(RP2str).SetTagVal(['chanEnable' RP2out],1);
	pause(.2);
elseif strcmp(cfg.sRP2_1circuit,'sphere_RP2_GWN')
	disp('Sounds are generated by RP2');
else
	disp('Sounds are not generated');
end

%% 
cfg.(RP2str).SetTagVal('delaySND1',snd.ondelay);  %
cfg.(RP2str).SetTagVal('delaySND2',snd.ondelay);  %
cfg.(RP2str).SetTagVal('soundDur1',dur);  % default?
cfg.(RP2str).SetTagVal('soundDur2',dur);  % default?


RP2out = cfg.mux2rp2(sndsetup(2)); 
	if strcmp(RP2str,'RP2_1')
		RP2out = 2; % in missing pa5 case, only output 2 of RP2.1 is used. 
	end
if strcmp(RP2str,'RP2_1')
	disp('PA5s 1 only; RP2.1');
	cfg.RA16_1.SetTagVal('rp2Enable1',1);
	cfg.RA16_1.SetTagVal('eventRP21',snd.onevent+1);
	if RP2out==1
			disp('PA5 1');
		cfg.PA5_1.SetAtten(atten);
	elseif RP2out==2
			disp('PA5 1');
		cfg.PA5_1.SetAtten(atten);
	end
elseif strcmp(RP2str,'RP2_2')
	disp('PA5s 2 and 3; RP2.2');
	cfg.RA16_1.SetTagVal('rp2Enable2',1);
	cfg.RA16_1.SetTagVal('eventRP22',snd.onevent+1);
	if RP2out==1
			disp('PA5 2');
		cfg.PA5_2.SetAtten(atten);
	elseif RP2out==2
			disp('PA5 3');
		cfg.PA5_3.SetAtten(atten);
	end
end

MUX(cfg.(RP2str),sndsetup(2),sndsetup(3));

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

% Stimulus-response plot+linear regression: azimuth
axes(handles.axes_regress_azimuth); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
unityline;
box off;
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('Target (deg)');
ylabel('Response (deg)');
title('Azimuth');

% Stimulus-response plot+linear regression: elevation
axes(handles.axes_regress_elevation); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
unityline;
box off
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('Target (deg)');
ylabel('Response (deg)');
title('Elevation');

axes(handles.axes_position); %#ok<*NASGU>
cla
hold on
box off
set(gca,'YTick',-90:30:90,'TickDir','out');
xlabel('Time (s)');
ylabel('Position (deg)');
ylim([-100 +100]);

axes(handles.axes_velocity); %#ok<*NASGU>
cla
hold on
box off
set(gca,'YTick',0:50:200,'TickDir','out');
xlabel('Time (s)');
ylabel('Velocity (deg)');
ylim([0 200]);

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

axes(handles.axes_xy); %#ok<*NASGU>
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
axes(handles.axes_xy);
handles.cfg.hcurdat(1) = plot(az,el,'k-','lineWidth',2);

axes(handles.axes_position);
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
disp('Waiting for RA16 button press/sound/led/acquisition');
% pause(.1); % short break
while ~cfg.RA16_1.GetTagVal('Wait')
	% 	disp('waiting')
	% do nothing
end

%% Sound play
% TODO: correct waiting/state machine
if strcmp(cfg.sRP2_1circuit,'RP2_WAV')
	disp('Waiting for RP2 sound');
	selsnd = strcmpi({stim.modality},'sound');
	if any(selsnd)
		snd		= stim(selsnd);
		nsnd	= numel(snd);
		for sndIdx = 1:nsnd
			sndsetup	= cfg.lookup(snd(sndIdx).Z+1,2:4);
			RP2out		= num2str(cfg.mux2rp2(sndsetup(2)));
			switch sndsetup(1)
				case 1
					cfg.RP2_1.GetTagVal(['bufferPos' RP2out])
					while cfg.RP2_1.GetTagVal(['bufferPos' RP2out])
						% do nothing
					end
				case 2
					cfg.RP2_2.GetTagVal(['bufferPos' RP2out])
					while cfg.RP2_2.GetTagVal(['bufferPos' RP2out])
						% do nothing
					end
			end
		end
	end
	
end

%% Data acquisition
disp('Waiting for data acquisition');
% pause(.1); % short break
while cfg.RA16_1.GetTagVal('Active')
	% do nothing
end

function data = trialSave(cfg,trial,data)
% DATA = SAVETRIAL(RA16_1,CFG,LED,DATA)
%
% Save data


%% Read RP2 DATA
% RP2_data	= NaN(cfg.nsamples,2);
% t = tic;
% for ii		= 1:cfg.nchan
% 	RP2_data(:,ii) = cfg.RP2_2.ReadTagV(cfg.recdataidx{ii},0,cfg.maxSamples+1000)';
% end
% rpdur			= toc(t);
% data(1).raw		= RA16_data;
% data(1).trialnr	= cfg.trial;
% data(1).RAsavedur	= radur;

%% Read DATA
RA16_data	= NaN(cfg.nsamples,cfg.nchan);
t = tic;
for ii		= 1:cfg.nchan
	RA16_data(:,ii) = cfg.RA16_1.ReadTagV(cfg.dataidx{ii},0,cfg.nsamples)';
end
radur			= toc(t);
data(1).raw		= RA16_data;
data(1).trialnr	= cfg.trial;
data(1).RAsavedur	= radur;

%% Read RA16 Event buffer
eventN					= cfg.RA16_1.GetTagVal('eventSize');
RA16_event				= cfg.RA16_1.ReadTagV('eventData',0,eventN)';
data(1).event				= 10^45*RA16_event/1.4; % to convert to scalars


%% Save trial data
% Because an experiment might be stopped at any time, we need to record
% data trial by trial in order to prevent data loss.
trialsingle		= trial(cfg.trial);
[~,fname,~]		= fileparts(cfg.fname); % remove extension
fname			= [fname '-' num2str(cfg.trial,'%04u')];
% sndname = fname;
fname			= fcheckext(fname,'sphere');
fname			= fullfile([cfg.dname filesep 'trial' filesep],fname);

save(fname,'data','cfg','trialsingle');

% sndname			= fcheckext(sndname,'rec');
% sndname			= fullfile([cfg.dname filesep 'trial' filesep],sndname);
% save(fname,'RP2_data');

function trialClean(stim,cfg)

%% Remove ledhandle if it exists
nstim = numel(stim);
for stmIdx = 1:nstim
	if isfield(stim(stmIdx),'ledhandle');
		if ~isempty(stim(stmIdx).ledhandle)
			stim(stmIdx).ledhandle.delete; 	% delete(leds)/switch off light;
		end
	end
end

%% Turn sounds off
% by switching off bits on inactive PM2Relay muliplexers
for muxIdx = 1:4
	MUX(cfg.RP2_1,muxIdx,0)
	MUX(cfg.RP2_2,muxIdx,0)
end

% and by setting buffersize to 0
for rpIdx = 1:2
	RPstr = ['RP2_' num2str(rpIdx)];
	cfg.RA16_1.SetTagVal(['rp2Enable' num2str(rpIdx)],0);
	for chanIdx = 1:2
		str		= ['bufferSize' num2str(chanIdx)];
		cfg.(RPstr).SetTagVal(str,0);
		cfg.(RPstr).SetTagVal(['chanEnable'  num2str(chanIdx)],0);
	end
end

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


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
