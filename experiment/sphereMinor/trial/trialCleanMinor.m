function trialCleanMinor(stim,cfg)

%% Remove ledhandle if it exists
nstim = numel(stim);
for stmIdx = 1:nstim
	if isfield(stim(stmIdx),'ledhandle')
		if ~isempty(stim(stmIdx).ledhandle)
			stim(stmIdx).ledhandle.delete; 	% delete(leds)/switch off light;
		end
	end
end

%% Turn sounds off
% by switching off bits on inactive PM2Relay muliplexers
for muxIdx = 1:2
	MUX(cfg.RZ6_1,muxIdx,0)
end

% % and by setting buffersize to 0
% for rpIdx = 1:2
% 	RPstr = ['RP2_' num2str(rpIdx)];
% 	cfg.RA16_1.SetTagVal(['rp2Enable' num2str(rpIdx)],0);
% 	for chanIdx = 1:2
% 		str		= ['bufferSize' num2str(chanIdx)];
% 		cfg.(RPstr).SetTagVal(str,0);
% 		cfg.(RPstr).SetTagVal(['chanEnable'  num2str(chanIdx)],0);
% 	end
% end

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
