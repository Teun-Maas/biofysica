function trialRunMinor(cfg,stim)
% RUNTRIAL(ZBUS,RA16)

%% Run zBus
cfg.zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
cfg.zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

%% Button Press
disp('Waiting for RZ6 button press/sound/led/acquisition');

t = tic;
while ~cfg.RZ6_1.GetTagVal('Wait')
% 	disp('waiting')
	% do nothing
end % sound is played after this loop exits.


%% Data acquisition

disp('Waiting for data acquisition');
pause(0.1)
while cfg.RZ6_1.GetTagVal('Active')
	% do nothing
	cfg.RZ6_1.GetTagVal('BufferPos_1')
end

duracq = toc(t);
disp(['Time elapsed: ' num2str(duracq) ' seconds']);

