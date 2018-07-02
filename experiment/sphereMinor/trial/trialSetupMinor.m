function [stim,cfg] = trialSetupMinor(cfg,stim)
% HLED = RUNTRIAL(RA16)
%
% Set up experimental parameters

%% Set TDT parameters
selled		= strcmpi({stim.modality},'LED') |  strcmpi({stim.modality},'SKY') |  strcmpi({stim.modality},'LASER') | strcmpi({stim.modality},'IRLED') ;
selacq		= strcmpi({stim.modality},'data acquisition');
% seltrg = strcmpi({stim.modality},'trigger');
selsnd		= strcmpi({stim.modality},'sound');

% selsndacq	= strcmpi({stim.modality},'sound acquisition');


%% LED
if any (selled)
	led		= stim(selled);
	nled	= numel(led);
	s		= ledpattern(nled*2); % Contains both ON and OFF events
	% Add timing cues
	for ledIdx = 1:nled
		strled1 = ['eventLED' num2str(2*ledIdx-1)];
		strled2 = ['eventLED' num2str(2*ledIdx)];
		cfg.RZ6_1.SetTagVal(strled1,led(ledIdx).onevent+1);
		cfg.RZ6_1.SetTagVal(strled2,led(ledIdx).offevent+1);
		
		strled1 = ['delayLED' num2str(2*ledIdx-1)];
		strled2 = ['delayLED' num2str(2*ledIdx)];
		cfg.RZ6_1.SetTagVal(strled1,led(ledIdx).ondelay);
		cfg.RZ6_1.SetTagVal(strled2,led(ledIdx).offdelay);
	end
	
	
	% Set colour and intensity
	if isfield(led,'colour')
		col = led(ledIdx).colour;
	else
		col = 1; % default to green
	end
	dum		= 0;
	for ii	= 1:2
		dum = dum+1;
		if ii==1
			s(dum).set(led(ledIdx).Z,cfg.ledcolours{col},1);
			s(dum).intensity(cfg.ledcolours{col},led(ledIdx).intensity); % sphere range 1-50
		else	% reset channels
			s(dum).set(led(ledIdx).Z,cfg.ledcolours{col},0);
		end
	end

% Write to stimulus
% use ledcontroller_pi class for raspberry_pi setup, as in the sphereMinor lab
stim(selled).ledhandle = ledcontroller_pi('dcn-led00','dcn-led01');
stim(selled).ledhandle.write(s);	

end

%% Acquisition
if any(selacq)
	acq	= stim(selacq);
	cfg.RZ6_1.SetTagVal('eventAcq',acq.onevent+1);
	cfg.RZ6_1.SetTagVal('delayAcq',acq.ondelay); % [ms]
	cfg.RZ6_1.SetTagVal('acqSamples',cfg.nsamples); % amount of DA samples
end

%% Sound
if any(selsnd)
	snd		= stim(selsnd);
	nsnd	= numel(snd);
	for sndIdx = 1:nsnd
		setSoundMinor(snd(sndIdx),cfg,'RZ6_1');
	end
end

if ~exist('maxSamples','var')
	maxSamples = 0;
end
cfg.maxSamples = maxSamples;


	
%% Wait for?
% This needs some tweaking
% search for latest event with longest offset
% which should also include sampling period and sound although this does not have an
% offevent
e				= [stim.offevent];
d				= [stim.offdelay];
mxevent			= max(e);
sel				= ismember(e,mxevent);
mxdelay			= max([d(sel) ceil(1000*cfg.nsamples./cfg.RZ6Fs) ]);

%%
cfg.RZ6_1.SetTagVal('eventWait',mxevent+1);
cfg.RZ6_1.SetTagVal('delayWait',mxdelay);


