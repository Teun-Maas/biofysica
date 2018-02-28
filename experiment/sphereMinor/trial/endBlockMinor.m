function endBlock(cfg)


snd.X			= 0;
snd.Y			= 0;
snd.channel		= [];
snd.detect		= [];
snd.event		= [];
snd.intensity	= 60;
snd.modality	= 'sound';
snd.offdelay	= 0;
snd.offevent	= 0;
snd.ondelay		= 0;
snd.onevent		= 0;
snd.matfile		= 'rehandel.mat';
snd.azimuth		= 3.5084e-15;
snd.elevation	= -3.5084e-15;
snd.Z			= 10;
snd.ledhandle	= [];

% load handel;
% y = resample(y,48828,Fs);
% snd = 
fname		= fullfile(cfg.snddir,snd.matfile);
if ~exist(fname,'file')
	disp('You do not have Handel installed.');
	disp('1000 monkeys are doing your work for you, and are now writing Handel.');
	load handel;
	snd = resample(y,48828,Fs);
	save(fname,'snd');
end
% setSound(snd,cfg,'RP2_1');
stim = trialSetup(cfg,snd);

trialRun(cfg,snd);
