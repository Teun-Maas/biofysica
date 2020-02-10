function pa_createmultiplesounds_tmp

%%
close all hidden
clear all hidden

% Periods of the sounds: 15   136   258   379   500 ms
% snd0xx = 15
% snd1xx = 136
% snd2xx = 258
% snd3xx = 379
% snd5xx = 500
% Dynamic onset: 500:100:1500 % so 11 IDs
% sndxx1 = 500
% sndxx2 = 600
% etc
% sndx11 = 1500 ms
% This will be made in pa_gensnd_avam
% snd000 = static only

onsetid = 1:11;
onset = 500:100:1500;
nonset = numel(onset);

periodid = [0:3 5];
period = linspace(15,500,5); % modulation period
nperiod = numel(period);
% vel,dens,md,durrip,durstat
for ii = 1:nperiod
	
	for jj = 1:nonset
		close all
		
		%% Default
		Fs		= 48828.125; % Freq (Hz)
		Nramp	= round(5/1000*Fs); % 5 ms
		Nlevel	= round(20/1000*Fs); % 5 ms
		eqfile = 'E:\MATLAB\PANDA\Donders\Setups\supplementary\visaton_equalizer.mat';
		
		col = bone(4);
		
		%% Step 1. Generate sound
		
		MF = round(1000./period(ii));
		durstat = onset(jj);
		durrip = 3000-durstat;
		snd = pa_genripple(MF,0,20,durrip,durstat,'display',true);
		
		snd = snd./max(abs(snd));
		
		subplot(211)
		plot(snd,'k-','Color',col(1,:),'LineWidth',4)
		hold on
		xlim([1 numel(snd)]);
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		box off
		title('Entire waveform');
		
		subplot(245)
		plot(snd(1:Nramp*3),'k-')
		hold on
		xlim([1 Nramp*3]);
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		title('Step 1. Signal generation');
		box off
		p		= audioplayer(snd,Fs);
		% 		playblocking(p);
		
		
		%% Step 2. Filter
		% pa_highpass
		% pa_lowpass
		% pa_equalizer
		snd = pa_equalizer(snd,eqfile);
		snd = snd./max(abs(snd));
		subplot(211)
		plot(snd,'k-','Color',col(2,:),'LineWidth',2)
		
		subplot(246)
		plot(snd(1:Nramp*3),'k-')
		xlim([1 Nramp*3]);
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		title('Step 2. Equalizer');
		box off
		% p		= audioplayer(snd,Fs);
		% playblocking(p);
		
		%% Step 3. Ramp
		% pa_ramp
		snd = pa_ramp(snd,Nramp);
		snd = snd./max(abs(snd));
		
		subplot(211)
		plot(snd,'k-','Color',col(3,:),'LineWidth',.5)
		
		subplot(247)
		plot(snd(1:Nramp*3),'k-')
		xlim([1 Nramp*3]);
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		title('Step 3. Ramp');
		box off
		% p		= audioplayer(snd,Fs);
		% playblocking(p);
		
		%% Step 4. 20 ms level ramp
		% pa_levelramp
		
		% snd = pa_fart_levelramp(snd,Nlevel);
		
		
		N = round(48828.125*(100+50)/1000);
		Nzero	= zeros(N,1); %This will create a vector with N zeros
		snd		= [Nzero;snd; Nzero]; %This pre- and ap-pends the ze
		
		snd = snd./max(abs(snd));
		
		subplot(211)
		plot(snd,'k-','Color',col(4,:),'LineWidth',.1)
		
		subplot(248)
		plot(snd(1:Nramp*3),'k-','LineWidth',1)
		xlim([1 Nramp*3]);
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		title('Step 4. Zero-prepend');
		box off
		% p		= audioplayer(snd,Fs);
		% playblocking(p);
		
		
		%% Step 5. save
		pa_datadir;
		% snd = zeros(size(snd));
		sndid = periodid(ii)*100+onsetid(jj);
		pa_writewav(snd,['snd' num2str(sndid,'%03i')]); % bb
		% pa_writewav(snd,'snd002'); % lp
		% pa_writewav(snd,'snd003'); % hp
	end
end
% print('-depsc','-painter',mfilename);


%% snd 000
close all;

%% Default
Fs		= 48828.125; % Freq (Hz)
Nramp	= round(5/1000*Fs); % 5 ms
Nlevel	= round(20/1000*Fs); % 5 ms
eqfile = 'E:\MATLAB\PANDA\Donders\Setups\supplementary\visaton_equalizer.mat';

col = bone(4);

% Step 1. Generate sound
snd = pa_genripple(0,0,20,0,3000,'display',true);
snd = snd./max(abs(snd));

subplot(211)
plot(snd,'k-','Color',col(1,:),'LineWidth',4)
hold on
xlim([1 numel(snd)]);
xlabel('Time (samples)');
ylabel('Amplitude (au)');
box off
title('Entire waveform');

subplot(245)
plot(snd(1:Nramp*3),'k-')
hold on
xlim([1 Nramp*3]);
xlabel('Time (samples)');
ylabel('Amplitude (au)');
title('Step 1. Signal generation');
box off
p		= audioplayer(snd,Fs);
% 		playblocking(p);


% Step 2. Filter
% pa_highpass
% pa_lowpass
% pa_equalizer
snd = pa_equalizer(snd,eqfile);
snd = snd./max(abs(snd));
subplot(211)
plot(snd,'k-','Color',col(2,:),'LineWidth',2)

subplot(246)
plot(snd(1:Nramp*3),'k-')
xlim([1 Nramp*3]);
xlabel('Time (samples)');
ylabel('Amplitude (au)');
title('Step 2. Equalizer');
box off
% p		= audioplayer(snd,Fs);
% playblocking(p);

% Step 3. Ramp
% pa_ramp
snd = pa_ramp(snd,Nramp);
snd = snd./max(abs(snd));

subplot(211)
plot(snd,'k-','Color',col(3,:),'LineWidth',.5)

subplot(247)
plot(snd(1:Nramp*3),'k-')
xlim([1 Nramp*3]);
xlabel('Time (samples)');
ylabel('Amplitude (au)');
title('Step 3. Ramp');
box off
% p		= audioplayer(snd,Fs);
% playblocking(p);

% Step 4. 20 ms level ramp
% pa_levelramp

% snd = pa_fart_levelramp(snd,Nlevel);


N = round(48828.125*(100+50)/1000);
Nzero	= zeros(N,1); %This will create a vector with N zeros
snd		= [Nzero;snd; Nzero]; %This pre- and ap-pends the ze

snd = snd./max(abs(snd));

subplot(211)
plot(snd,'k-','Color',col(4,:),'LineWidth',.1)

subplot(248)
plot(snd(1:Nramp*3),'k-','LineWidth',1)
xlim([1 Nramp*3]);
xlabel('Time (samples)');
ylabel('Amplitude (au)');
title('Step 4. Zero-prepend');
box off
% p		= audioplayer(snd,Fs);
% playblocking(p);


% Step 5. save
pa_datadir;
% snd = zeros(size(snd));
pa_writewav(snd,'snd000'); % bb
% pa_writewav(snd,'snd002'); % lp
% pa_writewav(snd,'snd003'); % hp