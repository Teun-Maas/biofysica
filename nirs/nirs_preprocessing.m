function nirs_preprocessing(varargin)
% NIRS_PREPROCESSING
%
% Preprocesses all NIRS data

%% Initialization
d				= keyval('dir',varargin,'/Users/marcw/DATA/Roos Cartignij/NIRS sessie/NH/');
% d				= keyval('dir',varargin,'/Users/marcw/DATA/Roos Cartignij/NIRS sessie/CI/');
cd(d);

%% Do stuff
expdir			= dir;
expnames		= {expdir.name};
idx				= vectorstrfind(expnames,'.');
expnames(idx)	= [];
nexp			= numel(expnames);

for expIdx = 1:nexp
	expname = expnames{expIdx};
	cd(expname)
	oxy = dir('*.oxy3');
	oxy = oxy.name;
	xls = dir('*order.xlsx');
	xls = xls.name;
	preprocess([d expname],oxy,xls);
	cd(d)
end

function preprocess(d,oxy,xls)
%% Events / stimuli
[stim,ustimtxt,ustim, stimV,stimA,stimAV]		= getnirseventcode_AVCI(d,xls);

%% Raw data
% read the raw data from file as one long continuous segment without any additional filtering
cfg					= [];
cfg.dataset			= oxy;
data_raw			= ft_preprocessing(cfg);

%% Timing / synchronisation
disp('Sync');

data_trigger			= nirs_trigger(data_raw,'method','ciav','chan',99:100,'fsdown',250);

idx						= numel(data_trigger.trigger{1});
data_trigger.stim		= stim(1:idx);
data_trigger.stimA		= stimA(1:idx); % 0 1 2 3 all stimuli
data_trigger.stimV		= stimV(1:idx); % 0 1 2 3
data_trigger.stimAV		= stimAV(1:idx); % 0 1 
data_trigger.ustimtxt	= ustimtxt; % unique stimulus text
data_trigger.ustim		= ustim; % unique stimulus codes


%% remove on and offset
a							= [data_trigger.trigger{1}];
t							= data_trigger.time{1}; % s
selt						= t>a(1)/data_trigger.fsample-10 & t<a(end)/data_trigger.fsample+10; % within 10 s of on- or offset
data_clean					= data_trigger;
b							= data_clean.time{1};
data_clean.time{1}			= b(selt);
b							= data_clean.trial{1};
b							= b(:,selt);
data_clean.trial{1}			= b;
b							= data_clean.sync{1};
b							= b(selt);
data_clean.sync{1}			= b;

dt							= [0 diff(selt)];
rmsample					= find(dt>0);
data_clean.trigger{1}		= data_clean.trigger{1}-rmsample;
data_clean.triggeroff{1}	= data_clean.triggeroff{1}-rmsample;

%% Resample
disp('Downsample');

cfg						= [];
cfg.resamplefs			= 10;
cfg.demean				= 'yes';
data_down				= ft_resampledata(cfg, data_clean);
% sync, trigger, triggertime
% data
data_down.sync{1}		= transpose(resample(transpose(data_clean.sync{1}),cfg.resamplefs,data_clean.fsample));
data_down.trigger{1}	= round(data_clean.trigger{1}/data_clean.fsample*cfg.resamplefs);
data_down.triggeroff{1}	= round(data_clean.triggeroff{1}/data_clean.fsample*cfg.resamplefs);

%% Hemodynamic prediction
disp('Hemo');

trg			= data_down.trigger{1}; % trigger indices for onset of the standard
trgoff		= data_down.triggeroff{1}; % trigger indices for onset of the standard
nstim		= numel(data_down.ustim);

t			= data_down.time{1};
mx			= numel(t);
stimsync	= zeros(nstim,mx);
HDpred		= stimsync;

for ii = 1:nstim
	sel						= data_down.stim==data_down.ustim(ii);
	trig					= trg(sel);
	trigoff					= trgoff(sel);
	for stdidx = 1:length(trig)
		begsample		= trig(stdidx) ;
		endsample		= trigoff(stdidx);
		stimsync(ii,begsample:endsample) = 1;
	end
end
for ii = 1:nstim
	X				= stimsync(ii,:);
	Y				= nirs_hdrfunction(1,X);
	HDpred(ii,:)	= Y;
end
X					= data_down.trial{1}(97,:)-data_down.trial{1}(100,:);
Y					= nirs_hdrfunction(1,X);
HDpred(nstim+1,:)	= Y;
data_down.HDR		= HDpred;

%% Scalp coupling index
disp('SCI');

data_sci				= nirs_sci(data_down);
% cfg = [];
% cfg.threshold = 0.7;
% data_sci				= ft_nirs_scalpcouplingindex(cfg,data_down);


%% remove heartbeat
disp('Heartbeat');

data_nophys				= nirs_removephysnoise(data_sci);

%% Bandpass
disp('Bandpass');

cfg						= [];
cfg.bpfilter			= 'yes';
cfg.bpfreq				= [0.05 0.4];
data_band				= ft_preprocessing(cfg,data_nophys);
data_nophys.trial{1}	= data_band.trial{1};

%% Transform
disp('Transform');

cfg						= [];
cfg.target				= {'O2Hb', 'HHb'};
cfg.channel				= 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all?, all NIRS channels are selected then
data_trans				= ft_transform_ODs(cfg, data_nophys);

fname = fcheckext([d filesep 'trans' oxy],'mat');
save(fname,'data_trans');

%% Reference channel subtraction
disp('RCS');
data_rcs				= nirs_rcs(data_trans); %#ok<*NASGU>

%% Motion artefacts
% disp('Motion Artefact');
% A						= data_rcs.trial{1};
% B						= movart2clean(A');
% B						= B';
% data_ma					= data_rcs;
% data_ma.trial{1}		= B;

%% Polynomial detrend
disp('Detrend');
data_detrend			= nirs_detrend(data_rcs);

%% Save
data					= data_detrend;
fname					= fcheckext([d filesep 'detrend' oxy],'mat');
save(fname,'data');

