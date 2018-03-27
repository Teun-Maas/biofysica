close all hidden
clearvars

%% Active X Control/Objects
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
RZ6_1circuit = which('sphereMinor_RZ6_soundrecording.rcx');

[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RZ6_1, err(2)]	= RZ6(1,RZ6_1circuit); % Real-time acquisition


%% TDT status
RZ6_1Status	= RZ6_1.GetStatus;

%% Run zBus
zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

pause(3)
%% Configuration

dur = 3; % s
nsamples = round(3*244.14)
% nsamples = 1000
RZ6_data		= NaN(nsamples,3);
tic;
for ii			= 1:3
	str = ['Data_' num2str(ii)]
	d = RZ6_1.ReadTagV(str,0,nsamples)';
	RZ6_data(:,ii) = d;

	toc
end
%%

t = (1:nsamples)/244.14;
figure(2)
subplot(311)
plot(t,RZ6_data(:,1));

subplot(312)
plot(t,RZ6_data(:,2));

subplot(313)
plot(t,RZ6_data(:,3));

%%
dur = 3; % s
nsamples = round(3*48828.125)
RZ6_mic		= NaN(nsamples,3);
for ii = 1:2
	str = ['Recording_' num2str(ii)];
	d = RZ6_1.ReadTagV(str,0,nsamples)';
	RZ6_mic(:,ii) = d;
	
	toc
end
%%

%%
t = (1:nsamples)/48828.125;
figure(3)
subplot(211)
plot(t,RZ6_mic(:,1));


%%
p = rms(RZ6_mic(:,1))
20*log10(p)
pref = 10^(94/20);
scaling = pref/p;
20*log10(scaling*p)

%%
figure(4)
clf
x = RZ6_mic(:,1);
[f,a] = getpower(x,48828.125,'display',true);

% [F,A,PH] = PA_GETPOWER(X,FS)
%%
Fs = 48828.125;
p = audioplayer(x,Fs);
play(p)
%%
% save ACOpacific_calibration RZ6_mic
save etymotic_speech RZ6_mic