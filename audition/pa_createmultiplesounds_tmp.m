function pa_createmultiplesounds_tmp

%%
close all hidden
clear all hidden

sndid = [100:199;200:299;300:399];
[ntyp,nrep] = size(sndid);

for ii = 1:ntyp
	for jj = 1:nrep
		close all
		
		%% Default
		Fs		= 48828.125; % Freq (Hz)
		Nramp	= round(5/1000*Fs); % 5 ms
		Nlevel	= round(20/1000*Fs); % 5 ms
		eqfile = 'E:\MATLAB\PANDA\Donders\Setups\supplementary\visaton_equalizer.mat';
		
		col = bone(4);
		
		%% Step 1. Generate sound
		% pa_gensweep - Schroeder sweep
		% pa_gengwn - Gaussian white noise from random time trace
		% pa_gengwnflat - Gaussian white noise from random phase spectrum
		% pa_gentone - Pure tone
		% pa_genripple - generate ripple
		
		snd = pa_gengwn(0.15);
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
		% playblocking(p);
		
		
		%% Step 2. Filter
		% pa_highpass
		% pa_lowpass
		% pa_equalizer
		snd = pa_equalizer(snd,eqfile);
		snd = snd./max(abs(snd));
		snd = pa_bandpass(snd);
		if ii == 2
			snd = pa_lowpass(snd);
		elseif ii == 3
			snd = pa_highpass(snd);
		end
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
		pa_writewav(snd,['snd' num2str(sndid(ii,jj),'%03i')]); % bb
		% pa_writewav(snd,'snd002'); % lp
		% pa_writewav(snd,'snd003'); % hp
	end
end
% print('-depsc','-painter',mfilename);
