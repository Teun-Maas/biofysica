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
% warning('spherePrime now checks each trial whether data should be acquired'); % this warning should be removed
% t = {trial(cfg.trial).stim.modality};
% if any(strcmp('data acquisition',t)) % check whether data should be acquired
	RA16_data		= NaN(cfg.nsamples,cfg.nchan);
	t				= tic;
	for ii			= 1:cfg.nchan
		RA16_data(:,ii) = cfg.RA16_1.ReadTagV(cfg.dataidx{ii},0,cfg.nsamples)';
	end
radur				= toc(t);
% else
% 		RA16_data		= [];
% 		radur = 0;
% end
data(1).raw			= RA16_data;
data(1).trialnr		= cfg.trial;
data(1).RAsavedur	= radur;

%% Read RA16 Event buffer
eventN					= cfg.RA16_1.GetTagVal('eventSize');
RA16_event				= cfg.RA16_1.ReadTagV('eventData',0,eventN)';
data(1).event			= 10^45*RA16_event/1.4; % to convert to scalars


%% Save trial data
% Because an experiment might be stopped at any time, we need to record
% data trial by trial in order to prevent data loss.
trialsingle		= trial(cfg.trial);
[~,fname,~]		= fileparts(cfg.fname); % remove extension
fname			= [fname '-' num2str(cfg.trial,'%04u')];
% sndname = fname;
fname			= fcheckext(fname,'sphere');
fname			= fullfile([cfg.dname filesep 'trial' filesep],fname);

%% remove graphics
if isfield(cfg,'hcurtar') 
	cfg				= rmfield(cfg,'hcurtar');
elseif isfield(cfg,'hcurdat') 
	cfg				= rmfield(cfg,'hcurdat');
else
	warning('No graphic handles present to be removed');
end

save(fname,'data','cfg','trialsingle');
