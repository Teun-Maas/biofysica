function pa_tdt_playsound(snd,RP2_1,zBus,com,level,speaker,stimSnd1,cmdSpeaker)% PA_TDT_PLAYSOUND%% Play sound in FART%% Dick Heeren / Tom / Denise% 2012 Modified: Marc van Wanrooij%% InitializationmaxSamples	= length(snd);stmSnd1		= 10;%% Write dataRP2_1.WriteTagV('WavData', 0, snd(1:maxSamples));%% Speaker onRP2_1.SetTagVal('Level',level); % Set speaker levelarg	= sprintf('%d;%d;1',speaker,stimSnd1); % micro argument, 1 - speaker onmicro_cmd(com,cmdSpeaker,arg); % give command to microcontroller%% Play soundzBus.zBusTrigA(0,0,0); % start sound% str		= sprintf('%d;%d;%d',2,2,1);% micro_cmd(com,100,str);busy = RP2_1.GetTagVal('Play'); % check sound playingwhile busy    busy = RP2_1.GetTagVal('Play');end % end sound%%	Stop sound playingRP2_1.SetTagVal('Level',0); % set level to 0arg	= sprintf('%d;%d;0',speaker,stimSnd1); % 0 - speaker offmicro_cmd(com,cmdSpeaker,arg); % give command to microcontroller% str		= sprintf('%d;%d;%d',2,2,0);% micro_cmd(com,100,str);