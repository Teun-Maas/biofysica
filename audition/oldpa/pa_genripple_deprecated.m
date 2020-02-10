function [snd,Fs] = pa_genripple(vel,dens,md,durrip,durstat,varargin)
% [SND,FS] = PA_GENRIPPLE(VEL,DENS,MOD,DURDYN,DURSTAT)
%
% Generate a ripple stimulus with velocity (amplitude-modulation) VEL (Hz),
% density (frequency-modulation) DENS (cyc/oct), and a modulation depth MOD
% (0-1). Duration of the ripple stimulus is DURSTAT+DURRIP (ms), with the first
% DURSTAT ms no modulation occurring.
%
% These stimuli are parametrized naturalistic, speech-like sounds, whose
% envelope changes in time and frequency. Useful to determine
% spectro-temporal receptive fields.  Many scientists use speech as
% stimuli (in neuroimaging and psychofysical experiments), but as they are
% not parametrized, they are basically just using random stimulation (with
% random sentences).  Moving ripples are a complete set of orthonormal
% basis functions for the spectrogram.
%
% See also PA_GENGWN, PA_WRITEWAV

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com
%
% Acknowledgements:
% Original script from Huib Versnel and Rob van der Willigen
% Original fft-script from NLS tools (Powen Ru)

%% Initialization
if nargin<1
	vel = 0; % (Hz)
end
if nargin<2
	dens = 1; % (cyc/oct)
end
if nargin<3
	md = 100; % Percentage (0-100%)
end
if nargin<4
	durrip = 1000; %msec
end
if nargin<5
	durstat = 500; %msec
end

%% Optional arguments
disp         = pa_keyval('display',varargin);
if isempty(disp)
	disp	= 1;
end
sv         = pa_keyval('save',varargin);
if isempty(sv)
	sv	= 'n';
end
plee         = pa_keyval('play',varargin);
if isempty(plee)
	plee	= 'y';
end
Fs         = pa_keyval('freq',varargin);
if isempty(Fs)
	Fs			= 48828.125; % Freq (Hz)
end
meth         = pa_keyval('method',varargin);
if isempty(meth)
	meth			= 'time';
	% other meth = 'fft' works slightly faster
end
% tic

md			= md/100; % Gain (0-1)
nRip        = round((durrip/1000)*Fs); % # Samples for Rippled Noise
nStat       = round((durstat/1000)*Fs); % # Samples for Static Noise
nTime		= nRip + nStat; % Total # of samples
time		= ((1:nTime)-1)/Fs; % Time (sec)

%% According to Depireux et al. (2001)
nFreq		= 128;
FreqNr		= 0:1:nFreq-1;
F0			= 250;
Freq		= F0 * 2.^(FreqNr/20);
Oct			= FreqNr/20;                   % octaves above the ground frequency
Phi			= pi - 2*pi*rand(1,nFreq); % random phase
Phi(1)		= pi/2; % set first to 0.5*pi


%% Generating the ripple
% Create amplitude modulations completely dynamic without static
T			= 2*pi*vel*time;
F			= 2*pi*dens*Oct;
[max(T) max(time) max(F) max(Oct) Fs/1000 nFreq durrip/1000 vel dens md]

[T,F]		= meshgrid(T,F);
A			= 1+md*sin(T+F);
A			= A';
whos A T F time Oct

% Modulate carrier, with static and dynamic part
snd		= 0;
for ii = 1:nFreq
	stat	= 1.*sin(2*pi* Freq(ii) .* time(1:nStat) + Phi(ii));
	rip		= A(1:nRip,ii)'.*sin(2*pi* Freq(ii) .* time(1:nRip) + Phi(ii));
	carr	= [stat rip];
	snd		= snd+carr;
end


%% Normalize ripple power to static power
rms_stat	= norm(snd(1:nStat)/sqrt(nStat));
rms_rip		= norm(snd(nStat+1:nStat+nRip)/sqrt(nRip));
ratio		= rms_stat/rms_rip;
% ratio		= 1;
snd			= [snd(1:nStat) ratio*snd(nStat+1:nStat+nRip)];

%% Normalization (so amplitude does not exceed 1)
snd			= snd/55; % in 3 sets of 500 trials, mx was 3 x 44+-1
% toc

%% Graphics
if disp
	plotspec(snd,Fs,durstat,nStat,A,Freq);
end

%% Play
if strcmpi(plee,'y');
	snd = pa_envelope(snd',round(10*Fs/1000));
	p	= audioplayer(snd,Fs);
	playblocking(p);
end

%% Save
if strcmpi(sv,'y');
	wavfname = ['V' num2str(vel) 'D' num2str(dens) '.wav'];
	pa_writewav(snd,wavfname);
end

function plotspec(snd,Fs,durstat,nStat,A,Freq)
close all;
t = (1:length(snd))/Fs*1000;
subplot(221)
plot(t,snd,'k-')
ylabel('Amplitude (au)');
ylabel('Time (ms)');
xlim([min(t) max(t)]);
hold on
plot(t+durstat,A(:,1)/2,'r-','LineWidth',2)

subplot(224)
pa_getpower(snd(nStat+1:end),Fs,'orientation','y');
ylim([min(Freq) max(Freq)])
ax = axis;
xlim(0.6*ax([1 2]));

subplot(223)
nsamples	= length(snd);
t			= nsamples/Fs*1000;
dt			= 12.5;
nseg		= t/dt;
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 1000;
spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
% colorbar
cax = caxis;
caxis([0.7*cax(1) 1.1*cax(2)])
ylim([min(Freq) max(Freq)])
set(gca,'YTick',[0.5 1 2 4 8 16]*1000);
set(gca,'YTickLabel',[0.5 1 2 4 8 16]);
set(gca,'YScale','log');
drawnow

function snd = genstat(durstat,Fs)
nTime       = round( (durstat/1000)*Fs ); % # Samples for Static Noise
time		= ((1:nTime)-1)/Fs; % Time (sec)

%% According to Depireux et al. (2001)
nFreq		= 128;
FreqNr		= 0:1:nFreq-1;
F0			= 250;
Freq		= F0 * 2.^(FreqNr/20);
Phi			= pi - 2*pi*rand(1,nFreq); % random phase
Phi(1)		= pi/2; % set first to 0.5*pi

%% Modulate carrier, with static and dynamic part
snd		= 0;
for ii = 1:nFreq
	stat	= 1.*sin(2*pi* Freq(ii) .* time + Phi(ii));
	snd		= snd+stat;
end