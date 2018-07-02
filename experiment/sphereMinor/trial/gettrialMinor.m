function handles = gettrialMinor(handles)


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

%% filter parameters
fname			= fullfile(cfg.expdir,cfg.expfname);
fname			= fcheckext(fname,'mat');
load(fname);
cfg.parameters = parameters;


%% CFG file
cfg						= readcfg(cfg); % read cfg cfile
cfg.acqdur				= cfg.humanv1.ADC(1).samples / cfg.humanv1.ADC(1).rate * 1000; % TODO: HumanV1/duration of data acquisition (ms)
cfg.nsamples			= round(cfg.acqdur/1000*cfg.RZ6Fs); % length data acquisition (samples)
cfg.nchan				= 3;


% TODO: code to verify experiment
trial					= sphereZMinor(trial,cfg);

%% handles
handles.trial	= trial;
handles.cfg		= cfg;
