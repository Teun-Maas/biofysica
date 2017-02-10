function maxSamples = setSound(snd,cfg,RP2str)
% SETSOUND(SND,CFG,RP2STR)
%
% Set all parameters for sound presentation

Z			= snd.Z;
atten		= max(cfg.maxsndlevel-snd.intensity,0);

dur			= snd.offdelay-snd.ondelay;
sndsetup	= cfg.lookup(Z+1,2:4);


%% WAV
if strcmp(cfg.sRP2_1circuit,'sphere_RP2_WAV')
	disp('Sounds are obtained from MAT file');
	fname		= fullfile(cfg.snddir,snd.matfile);
	s			= load(fname); % Read sound from mat file
	sig			= transpose(s.snd);
	maxSamples	= length(sig);
	maxSamples	= min([maxSamples,500000]); % maximum set in RP2 circuit
	dur			= round(1000*maxSamples/cfg.medusaFs);
	RP2out		= num2str(cfg.mux2rp2(sndsetup(2))); % from MUX to RP2 channel
	cfg.(RP2str).SetTagVal(['bufferSize' RP2out],maxSamples);
	cfg.(RP2str).WriteTagV(['wavData' RP2out], 0,sig(1:maxSamples));
	cfg.(RP2str).SetTagVal(['chanEnable' RP2out],1);
	pause(.2);
elseif strcmp(cfg.sRP2_1circuit,'sphere_RP2_GWN')
	disp('Sounds are generated by RP2');
else
	disp('Sounds are not generated');
end

%% 
cfg.(RP2str).SetTagVal('delaySND1',snd.ondelay);  %
cfg.(RP2str).SetTagVal('delaySND2',snd.ondelay);  %
cfg.(RP2str).SetTagVal('soundDur1',dur);  % default?
cfg.(RP2str).SetTagVal('soundDur2',dur);  % default?


RP2out = cfg.mux2rp2(sndsetup(2));
if strcmp(RP2str,'RP2_1')
	disp('PA5s 1 and 3; RP2.1');
	disp(['MUX = ' num2str(sndsetup(2))]);
	disp(['RP2 channel = ' num2str(RP2out)]);
	cfg.RA16_1.SetTagVal('rp2Enable1',1);
	cfg.RA16_1.SetTagVal('eventRP21',snd.onevent+1);
	if RP2out==1
		disp('PA5 1');
		cfg.PA5_1.SetAtten(atten);
	elseif RP2out==2
		disp('PA5 3');
		cfg.PA5_3.SetAtten(atten);
	end
elseif strcmp(RP2str,'RP2_2')
	disp('PA5s 2 and 4; RP2.2');
	disp(['MUX = ' num2str(sndsetup(2))]);
	disp(['RP2 channel = ' num2str(RP2out)]);
	cfg.RA16_1.SetTagVal('rp2Enable2',1);
	cfg.RA16_1.SetTagVal('eventRP22',snd.onevent+1);
	if RP2out==1
			disp('PA5 2');
		cfg.PA5_2.SetAtten(atten);
	elseif RP2out==2
			disp('PA5 4');
		cfg.PA5_4.SetAtten(atten);
	end
end

MUX(cfg.(RP2str),sndsetup(2),sndsetup(3));
