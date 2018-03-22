function trialRunMinor(cfg,stim)
% RUNTRIAL(ZBUS,RA16)

%% Run zBus
cfg.zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
cfg.zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

%% Button Press
disp('Waiting for RZ6 button press/sound/led/acquisition');
% pause(.1); % short break
while ~cfg.RZ6_1.GetTagVal('Wait')
	% 	disp('waiting')
	% do nothing
end

%% Sound play
% TODO: correct waiting/state machine
if strcmp(cfg.sRZ6_1circuit,'RP2_WAV')
	disp('Waiting for RZ6 sound');
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
while cfg.RZ6_1.GetTagVal('Active')
	% do nothing
end
