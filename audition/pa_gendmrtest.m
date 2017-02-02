function [snd,Fs] = pa_gendmr(vel,dens,md,durrip,durstat,varargin)
% [SND,FS] = PA_GENDNR(VEL,DENS,MOD,DURDYN,DURSTAT)
%
% Generate a Dynamic Moving Ripple
%
% SSDMR(t,Xk) = M/2 * sin[2 * pi * Omega(t) * Xk + Phi(t)]
% M = modulation depth (dB)
% Xk = log2(fk/f1) octave frequency axis relative to f1
% Phi(t) = integral 0-t Fm(tau)d(tau) controls time-varying temporal
% modulation rate Fm(t)
% Spectral [Omega(t)] and temporal [Fm(t)] parameters are independent,
% slowyly timevarying random processes with maximum rates of change 1.5 Hz
% for Fm and 3 Hz for Omega (maximum ranges in speech and vocalizations,
% Greenberg 1998)
% uniformly flat distributed amplitudes 0–4 cycles per octave for ? and
% ?350 to +350 Hz for Fm
%
clc
close all hidden;
clear all hidden;

%% Initialization
dur				= 3; % s
velocityrate	= 3; % temporal rate of change (Hz)
densityrate		= 6; % spectral rate of change (Hz)
ndensity		= round(dur*densityrate);
nvelocity		= round(dur*velocityrate);
Fs				= 44100; % sampling rate (Hz)
densityrange	= [0 4]; % (cycles/octave)
velocityrange	= [-64 64]; % (Hz)
densityrange	= [0 2]; % (cycles/octave)
velocityrange	= [-32 32]; % (Hz)

%% Create envelope
density		= randn(ndensity,1);
velocity	= randn(nvelocity,1);

subplot(221)
t = (0:length(density)-1)/densityrate;
plot(t,density,'ko-')
hold on

subplot(222)
t = (0:length(velocity)-1)/velocityrate;
plot(t,velocity,'ko-')
hold on
% plot((0:length(velocity)-1)/velocityrate,velocity,'ro-');

%% upsample
% resample
nd			= Fs/densityrate;
density		= resample(density,nd,1);
nd = numel(density);

nv			= Fs/velocityrate;
velocity		= resample(velocity,nv,1);
nv = numel(velocity);

N = min([nd nv]);

density = density(1:N);
velocity = velocity(1:N);

subplot(221)
t = (0:length(density)-1)/Fs;
plot(t,density,'r-')

subplot(222)
t = (0:length(velocity)-1)/Fs;
plot(t,velocity,'r-')

% interpolation
% interp1(x,y,xi)
subplot(223)
hist(density,50);

subplot(224)
hist(velocity,50);

%% Convert to uniform distribution
velocity	= erf(velocity);
density		= erf(density);

%% Set range
velocity	= velocity*range(velocityrange)/2+mean(velocityrange);
density		= density*range(densityrange)/2+mean(densityrange);

% density = zeros(size(density));
% velocity = repmat(10,size(velocity));

figure
subplot(223)
hist(density,50);

subplot(224)
hist(velocity,50);

%% Generating the ripple
% SSDMR(t,Xk) = M/2 * sin[2 * pi * Omega(t) * Xk + Phi(t)]
% M = modulation depth (dB)
% Xk = log2(fk/f1) octave frequency axis relative to f1
% Phi(t) = integral 0-t Fm(tau)d(tau) controls time-varying temporal
% modulation rate Fm(t)
M				= 45;
time			= (0:length(density)-1)/Fs;
time = time';
nFreq			= 2^7;
FreqNr			= (0:1:nFreq-1)/nFreq;
F0				= 500;
Fmax			= 16000;
Orange = pa_freq2bw(F0,Fmax);
Xk				= Orange*FreqNr;                   % octaves above the ground frequency

%% Try repmat instead of for-loop
fk		= pa_oct2bw(F0,Xk);
phi		= 2*pi*rand(size(fk));
velocity = 2*pi*cumsum(velocity/Fs);

snd = zeros(size(time));
for ii	= 1:nFreq
	Sdb		= 1+sin(2*pi*density.*Xk(ii)+velocity);
	% 	Slin	= 10.^(Sdb-M/2)/20;
	% 	Slin = Sdb./max(Sdb);
	Slin = Sdb;
	snd		= snd+Slin.*sin(2*pi*fk(ii)*time+phi(ii));
end
snd = snd/max(snd);
figure
ax(1) = subplot(211);
plot(time,snd)


density		= repmat(density,1,size(Xk,2));
Xk			= repmat(Xk,size(density,1),1);
velocity	= repmat(velocity,1,size(Xk,2));
time		= repmat(time,1,size(Xk,2));
Sdb		= M/2*sin(2*pi*density.*Xk+2*pi*velocity.*time);
ax(2) = subplot(212);
imagesc(time(:,1),Xk(:,1),Sdb')
drawnow

linkaxes(ax,'x');

figure;
imagesc(time(:,1),Xk(:,1),Sdb');
set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Frequency (Octaves)');
axis square;
xlim([0 0.6]);
% keyboard

%% Save graphics
pa_datadir;
print('-depsc2','-painter',mfilename);

%% Save sound
pa_writewav(snd,mfilename,0.99,Fs);
wavplay(snd,Fs)
