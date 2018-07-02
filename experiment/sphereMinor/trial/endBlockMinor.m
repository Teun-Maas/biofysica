function endBlockMinor(cfg)


stim.X			= 0;
stim.Y			= 0;
stim.channel		= [];
stim.detect		= [];
stim.event		= [];
stim.intensity	= 60;
stim.modality	= 'sound';
stim.offdelay	= 0;
stim.offevent	= 0;
stim.ondelay	= 0;
stim.onevent	= 0;
stim.matfile	= 'rehandel.mat';
stim.azimuth	= 3.5084e-15;
stim.elevation	= -3.5084e-15;
stim.Z			= 10;
stim.ledhandle	= [];
stim.parameters	= 99;

fname		= fullfile('C:\DATA\SND',stim.matfile);
if ~exist(fname)
	disp('You do not have Handel installed.');
	disp('1000 monkeys are doing your work for you, and are now writing Handel.');
	load handel;
	snd = resample(y,round(cfg.RZ6Fs),Fs); %#ok<NASGU>
	nsamples = numel(snd);
	dur = round( nsamples/cfg.RZ6Fs*1000);
	stim.duration = dur;
	save(fname,'snd');
end

disp('Hallelujah');
warning('cannot load into RZ6 circuit as it cannot handle .mat files');
% stim = trialSetupMinor(cfg,stim);
% trialRunMinor(cfg,stim);

%% Let's run some LEDs instead
import org.zeromq.ZMQ

n = 5;
s = ledpattern(n);

ir=50;
ig=ir;
for i=1:n
	if mod(i,2) == 0
		s(i).set(0:2:31,'r');
	else 
		s(i).set(1:2:31,'g');
	end
    s(i).intensity('r', ir);
    s(i).intensity('g', ig);
end

leds = ledcontroller_pi('dcn-led00','dcn-led01');

leds.write(s);
for i=1:n
    leds.trigger;
	pause(0.7);
end

%% Mop up
% Turn off the lights
t=ledpattern; % empty ledpattern
leds.write(t);
leds.trigger;
% delete objects
delete(leds);
delete(s);
delete(t);

