function data = trialSaveMinor(cfg,trial,data)
% DATA = SAVETRIAL(RA16_1,CFG,LED,DATA)
%
% Save data

%% Read DATA
% t = {trial(cfg.trial).stim.modality};
% if any(strcmp('data acquisition',t)) % check whether data should be acquired
RZ6_data		= NaN(cfg.nsamples,cfg.nchan);
t				= tic;

for	ii			= 1:cfg.nchan % only 3 data channels...
	RZ6_data(:,ii) = cfg.RZ6_1.ReadTagV(cfg.dataidx{ii},0,cfg.nsamples)';
end
dur				= toc(t);
data(1).raw			= RZ6_data;

%%
data(1).trialnr		= cfg.trial;
data(1).savedur		= dur;

stim = trial.stim;

selsnd		= strcmpi({stim.modality},'sound');
if any(selsnd) % Check if a sound was played
	stim		= stim(selsnd);
	stim		= stim(1); % stores only 1 sound
	nsamples = round(stim.duration/1000*cfg.RZ6Fs);
	snd = cfg.RZ6_1.ReadTagV('sound',0,nsamples)'; %#ok<*NASGU>
else % Defaults when no sound was played
	stim = nan;
	nsamples = 0;
	snd = nan;
end

%% Save trial data
% Because an experiment might be stopped at any time, we need to record
% data trial by trial in order to prevent data loss.
trialsingle		= trial(cfg.trial);
[~,fname,~]		= fileparts(cfg.fname); % remove extension
fname			= [fname '-' num2str(cfg.trial,'%04u')];
% sndname = fname;
fname			= fcheckext(fname,'sphere');
fname			= fullfile([cfg.dname filesep 'trial' filesep],fname);

%% remove graphics and objects
% from saved data, not from GUI handles
delhandles	= {'hcurtar','hcurdat','RZ6_1','zBus'};
ndelhandles = numel(delhandles);
for ii = 1:ndelhandles
	if isfield(cfg,delhandles{ii})
		cfg				= rmfield(cfg,delhandles{ii});
	end
end
save(fname,'data','cfg','trialsingle','snd');
