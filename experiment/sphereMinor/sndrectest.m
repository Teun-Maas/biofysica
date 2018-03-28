close all
clearvars

load('/Users/marcw/Dropbox/manuscripts/Divers/Patientlab/ACOpacific_calibration.mat');
Fs			= 48828.125;
x			= RZ6_mic(:,1);
nsamples	= length(x);
t			= (1:nsamples)/Fs;
p			= rms(x);
pref		= 10^(94/20);
scaling		= pref/p;
y			= scaling*x;


figure(1)
subplot(211)
plot(t,y);
figure(2)
clf
[f,a] = getpower(y,48828.125,'display',true,'nfft',1024);
ylim([-10 110]);
title(20*log10(rms(y)))

% p = audioplayer(x,Fs);
% play(p)



%%
load('/Users/marcw/Dropbox/manuscripts/Divers/Patientlab/ACOpacific_labnoise.mat');
Fs			= 48828.125;
x			= RZ6_mic(:,1);
nsamples	= length(x);
t			= (1:nsamples)/Fs;
y			= scaling*x;

figure(1)
subplot(212)
plot(t,y);

figure(2)

[f,a] = getpower(y,48828.125,'display',true,'nfft',1024,'Color','r');
ylim([-10 110]);
xlim([100 15000]);

legend('1000 Hz @ 94 dB','Background noise');
ylabel('Sound Pressure Level (dB), A-weighted');
%%
cd('/Users/marcw/Dropbox/manuscripts/Divers/Patientlab');
savegraph('background','eps');


%%
load('/Users/marcw/Dropbox/manuscripts/Divers/Patientlab/ACOpacific_speech.mat');
Fs			= 48828.125;
x			= RZ6_mic(:,1);
nsamples	= length(x);
t			= (1:nsamples)/Fs;
y			= scaling*x;


% figure(1)
% subplot(211)
% plot(t,y);
% figure(2)
% clf
% [f,a] = getpower(y,48828.125,'display',true,'nfft',1024);
% ylim([-10 110]);
% title(20*log10(rms(y)))

p = audioplayer(x,Fs);
play(p)
