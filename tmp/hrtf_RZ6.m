function hrtf_RZ6
%
%
%

close all hidden
% clearvars

%% Active X Control/Objects
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
RZ6_1circuit	= which('hrtf_RZ6.rcx');
zBus			= ZBUS(1); % zBus, number of racks
RZ6_1			= RZ6(1,RZ6_1circuit); % Real-time acquisition

%% Run zBus
zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).


%% Default recording (check in circuit to see whether this is true)
nsamples	= 2048*20;
%%
% dev 1 = B - 11-14 is -10:-10:-40 deg
% dev 2 = A - 11-16 is 10:10:60 deg

dev = [ones(1,4) zeros(1,7)]+1;
ch	= [14:-1:11 1 11:16];
nch = numel(ch);

mic		= NaN(nch,nsamples);
stim	= mic;
framp = mic;
for ii = 1:nch
	[mic(ii,:),stim(ii,:),framp(ii,:)]	= playsound(RZ6_1,zBus,nsamples,dev(ii),ch(ii));
		pause(0.1);
		zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).
		pause(0.1);
end

zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).


%%
t = (1:nsamples)/48828.125;
figure(3)
plot(t,mic);

%%
figure(4)
clf
x = mic(:);
[f,a] = getpower(x,48828.125,'display',true,'nfft',2^10);

% [F,A,PH] = PA_GETPOWER(X,FS)
%%
% Fs = 48828.125;
% p = audioplayer(x,Fs);
% play(p)
%%
save garcia_left_hrtf_6 mic stim framp
% save ACOpacific_labnoise RZ6_mic

% 1 = 60 dB amp
% 2 = 50 dB Amp
% 3 = head at +5 
% 4 = head at +40 
% 5 = head at -30 
% 6 = head at +15 

% save etymotic_speech RZ6_mic

function [rec,stim,framp] = playsound(RZ,zBus,nsamples,dev,ch)
%% Set speaker
MUX(RZ,dev,ch);
% pause(1)

%% Trigger Sound Generation
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

busy		= 1;

tic
pause(nsamples/48828.125+0.1);
while busy
	busy	= RZ.GetTagVal('Busy');
end
toc
%%
rec		= RZ.ReadTagV('Recording',0,nsamples)';
stim	= RZ.ReadTagV('Playing',0,nsamples)';
framp	= RZ.ReadTagV('FreqRamp',0,nsamples)';
