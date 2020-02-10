function pa_createmultiplesounds_test

close all hidden
sndid = [100 200 300 400 1 2 3 4 5 6 7 8 9 10]; % soundids, bb, hp, lp, tones kHz
for indx = sndid
	figure(indx)
	%% Default
	Fs		= 48828.125; % Freq (Hz)
	Nramp	= round(5/1000*Fs); % 5 ms
	Nlevel	= round(20/1000*Fs); % 5 ms
	dur		= 0.15; % s
	eqfile = 'E:\MATLAB\PANDA\Donders\Setups\supplementary\visaton_equalizer.mat';
	
	col = bone(4);
	
	%% Step 1. Generate sound
	% pa_gensweep - Schroeder sweep
	% pa_gengwn - Gaussian white noise from random time trace
	% pa_gengwnflat - Gaussian white noise from random phase spectrum
	% pa_gentone - Pure tone
	% pa_genripple - generate ripple
	
	if ismember(indx,[100 200 300])
		snd = pa_gengwn(dur);
	else
		snd = pa_gentone(indx*1000,dur*1000);
	end
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
% 	p		= audioplayer(snd,Fs);
	% playblocking(p);
	
	
	%% Step 2. Filter
	snd = pa_equalizer(snd,eqfile);
	snd = snd./max(abs(snd));

	if indx==300
		snd = pa_bandpass(snd);
		snd = pa_lowpass(snd);		
	elseif indx==400
		snd = pa_lowpass(snd);
	elseif indx==200
		snd = pa_highpass(snd);
	elseif indx==100
		snd = pa_bandpass(snd);
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
	snd = pa_fart_levelramp(snd,Nlevel);
	
	
	N		= round(48828.125*(100+50)/1000);
	Nzero	= zeros(N,1); %This will create a vector with N zeros
	snd		= [Nzero;snd; Nzero]; %#ok<AGROW> %This pre- and ap-pends the ze
	
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

	
	
	%% Step 5. save
	pa_datadir;
	cd('SND');
	pa_writewav(snd,['snd' num2str(indx,'%03i')]); % bb
	drawnow
end
