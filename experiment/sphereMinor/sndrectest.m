close all
clearvars

load('/Users/marcw/Dropbox/Work/Divers/Patientlab/ACOpacific_calibration.mat');
Fs			= 48828.125;
x			= RZ6_mic(:,1);
nsamples	= length(x);
t			= (1:nsamples)/Fs;
p			= rms(x);
pref		= 10^((94-0.3)/20);
scaling		= pref/p
y			= scaling*x;
% y			= filterA(y,Fs);

figure(1)
subplot(211)
plot(t,x);
figure(2)
clf
[f,a] = getpower(y,48828.125,'display',true,'nfft',2^12);
ylim([-10 110]);
title(20*log10(rms(y)));

w = aweight(f);
m = a+w;
		
p		= 10.^(a/20);
psum	= sqrt(sum(p.^2));
L		= 20*log10(psum);

p		= 10.^(m/20);
psum	= sqrt(sum(p.^2));
La		= 20*log10(psum);

levels = [20*log10(rms(y)) L La]

str = ['L_{rms} = ' num2str(round(levels(1)))];
text(1000+100,94-3,str,'FontSize',20);

str = ['L_{spec} = ' num2str(round(levels(2)))];
text(1000+100,94-6,str,'FontSize',20);

str = ['L_{spec,A} = ' num2str(round(levels(3)))];
text(1000+100,94-9,str,'FontSize',20);
% p = audioplayer(x,Fs);
% play(p)



%%
load('/Users/marcw/Dropbox/Work/Divers/Patientlab/ACOpacific_labnoise.mat');
Fs			= 48828.125;
x			= RZ6_mic(:,1);
nsamples	= length(x);
t			= (1:nsamples)/Fs;
y			= scaling*x;
% y = filterA(y,Fs);

figure(1)
subplot(212)
plot(t,x);

figure(2)

[f,a] = getpower(y,48828.125,'display',true,'nfft',2^12,'Color','r');
p		= 10.^(a/20);
psum	= sqrt(sum(p.^2));
L		= 20*log10(psum);

w = aweight(f);
m = a+w;

p		= 10.^(m/20);
psum	= sqrt(sum(p.^2));
La		= 20*log10(psum);

levels = [20*log10(rms(y)) L La];

str = ['L_{rms} = ' num2str(round(levels(1)))];
text(1000-100,30-3,str,'FontSize',20,'Color','r');

str = ['L_{spec} = ' num2str(round(levels(2)))];
text(1000-100,30-6,str,'FontSize',20,'Color','r');

str = ['L_{spec,A} = ' num2str(round(levels(3)))];
text(1000-100,30-9,str,'FontSize',20,'Color','r');

ylim([-10 110]);
xlim([100 15000]);
horline(94);
legend('1000 Hz @ 94 dB','Background noise');
ylabel('Sound Pressure Level (dB), A-weighted');
%%
% cd('/Users/marcw/Dropbox/manuscripts/Divers/Patientlab');
% savegraph('background','eps');


%%
load('/Users/marcw/Dropbox/Work/Divers/Patientlab/ACOpacific_speech.mat');
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
