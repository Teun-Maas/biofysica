function handles = getcfgdefaultsMinor(handles)
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
cfg = tdt_globalsMinor(cfg);		% TDT defaults, in cfg structure, includes circuit names for [RA16_1,RA16_2,RX6,RP2_1,RP2_2]

%% TDT and PLC defaults
cfg.nleds				= 8; % maximum number of LED configurations on PLC

cfg.calfile				= which('sphere.net');
if isempty(cfg.calfile)
	cfg.calfile				= which('defaultsphereMinor.net');
end

%% Convert experimental trial

cfg						= spherelookupMinor(cfg); % lookup structure

%% Data filter coefficients
% Do we use this
cfg.Fcutoff             = 80;       % Cutoff Frequency  (Hz)
cfg.Order                   = 50;

cfg.lpFilt = designfilt('lowpassfir', 'FilterOrder', cfg.Order, 'CutoffFrequency', ...
                    cfg.Fcutoff, 'SampleRate', cfg.RZ6Fs, 'Window', 'hamming');
				% check RZ6Fs
%% Led colours
cfg.ledcolours = {'g','r'}; % green and red

%% handles
handles.cfg		= cfg;
