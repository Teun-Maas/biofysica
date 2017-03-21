function stim = trialSetup(cfg,stim)
% HLED = RUNTRIAL(RA16)
%
% Set up experimental parameters

%% Set TDT parameters
selled		= strcmpi({stim.modality},'LED') |  strcmpi({stim.modality},'SKY') |  strcmpi({stim.modality},'LASER') | strcmpi({stim.modality},'IRLED') ;
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

	%%
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
		if isfield(led,'colour')
			col = led(ledIdx).colour;
		else
			col = 1;
		end
		for ii	= 1:2
			cnt = cnt+1;
			if ii==1
				s(cnt).set(led(ledIdx).Z,cfg.ledcolours{col},1);
				s(cnt).intensity(cfg.ledcolours{col},led(ledIdx).intensity); % hoop: range 0-255, sphere range 1-50
			else
				s(cnt).set(led(ledIdx).Z,cfg.ledcolours{col},0);
			end
		end
	end
	stim(find(selled,1)).ledhandle = ledcontroller;
	stim(find(selled,1)).ledhandle.write(s);
	%%
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
mxevent			= max(e);
sel				= ismember(e,mxevent);
mxdelay			= max([d(sel) ceil(1000*cfg.nsamples./cfg.medusaFs) ceil(1000*maxSamples/48828.125)]);

%%
cfg.RA16_1.SetTagVal('eventWait',mxevent+1);
cfg.RA16_1.SetTagVal('delayWait',mxdelay);
