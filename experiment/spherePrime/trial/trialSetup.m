function [stim,cfg] = trialSetup(cfg,stim)
% HLED = RUNTRIAL(RA16)
%
% Set up experimental parameters

% GW/20191030 Modifications for using LED patterns.

%% Set TDT parameters
selled		= strcmpi({stim.modality},'LED') |  strcmpi({stim.modality},'SKY') |  strcmpi({stim.modality},'LASER') | strcmpi({stim.modality},'IRLED') ;
selacq		= strcmpi({stim.modality},'data acquisition');
% seltrg = strcmpi({stim.modality},'trigger');
selsnd		= strcmpi({stim.modality},'sound');
% selsndacq	= strcmpi({stim.modality},'sound acquisition');


%% LED
if any(selled)
% 	nled	= numel(led);
% 	% 	nled = 2
% 	n		= nled*2; % LEDs need to be turned on and off
% 	s		= ledpattern(n);
    [patterns, pattern_times, pattern_events] = makeOrderedLedPatterns(stim,cfg);
    nEvents = length(pattern_times);

    nEventsmax = 8; %% RA16 circuit can handle a maximum of 8 events

    for ii = 1:nEventsmax
        % write events and times to TDT RA16
        % Set timing information on LEDs
        % Note that in RA16 circuit, event 1 = start of experiment
        str1 = ['eventLED' num2str(ii)];
        str2 = ['delayLED' num2str(ii)];
        if ii <= nEvents
            cfg.RA16_1.SetTagVal(str1, events(ii) + 1); %% event numbers in RA16 circuit are 1:8 instead of 0:7
            cfg.RA16_1.SetTagVal(str2, delays(ii) + 1); %% zero is not accepted by the RA16 circuit. All delay times are offset by 1 ms.
        else
            cfg.RA16_1.SetTagVal(str1, -1);             %% -1 is no event
            cfg.RA16_1.SetTagVal(str2, 1);              %% 1 is no delay; zero or -1 is not accepted by the RA16 circuit. 
        end
    end

    ledbox = ledcontroller;
    ledbox.write(patterns);
	stim(find(selled,1)).ledhandle = ledbox;
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
cfg.maxSamples = maxSamples;

%% Sound Acquisition
% if any(selsndacq)
% 	sndacq	= stim(selsndacq);
% 	cfg.RP2_1.SetTagVal('eventAcq',sndacq.onevent+1);
% 	cfg.RP2_1.SetTagVal('delayAcq',sndacq.ondelay);
% 	cfg.RP2_1.SetTagVal('acqSamples',cfg.nsamples); % amount of DA samples
% end

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
