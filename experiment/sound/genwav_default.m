function genwav_default
% GENWAV_DEFAULT
%
% Generate 300 wav files for the SPHERE setup:
% - 100 broadband Gaussian White Noises
% - 100 High-pass @ 3 kHz
% - 100 Low-Pass @ 1.5 kHz
%
% Wav files are stored according to convention:
% snd###.wav with ### 100-199 for BB, 200-299 for HP, 300-399 for LP
%
% Difference with HOOP: no 20 ms pre-zeros
%
% See also GENGWN, GENEXP_DEFAULTLOC

%% Initialization
close all hidden
clear hidden

%% Default
Fs		= 48828.125; % Default freq (Hz) of TDT RP2
Nramp	= round(5/1000*Fs); % # samples of cosine-squared on-offset ramp of 5 ms
dur		= 0.15; % sound duration = 0.15 s
Fc		= 20000; % maximum frequency
Fh		= 500; % high-pass @ 500 Hz
% freq = {'BB','HP','LP'}; % broadband, high-pass, low-pass
hp		= 3000; % high-pass (ILD and spectral cues)
lp		= 1500; % low-pass (ITD cues)
order	= 100; % arbitrary high filter order
for freqIdx = 1:3
	for sndIdx = 1:100
		snd		= gengwn(dur,order,Fc,Fs,Fh);
		snd		= equalizer(snd);
		switch freqIdx
			case 2
				snd = highpass(snd,'Fc',hp,'Fs',Fs,'order',order); % high-pass filter
			case 3
				snd = lowpass(snd,'Fc',lp,'Fs',Fs,'order',order); % low-pass filter
		end
		
		%% rms scaling
% 		P(freqIdx,sndIdx) = peak2rms(snd);
		snd		= 10*snd./(5*rms(snd)); % normalize power with 4*RMS(for HOOP, normalization was done on max between -1 and 1)
% 		M(freqIdx,sndIdx) = max(abs(snd));
		S(freqIdx,sndIdx) = sum(abs(snd)>10);
		sel = snd>10;
		snd(sel) = 10;
		sel = snd<-10;
		snd(sel) = -10;
		
		%% max scaling
		
		snd		= ramp(snd,Nramp); % on- and offset ramp
	% store in datadir
	filename = fullfile(pwd,['snd' num2str(sndIdx-1+100*freqIdx,'%03u') '.mat']); % convention: snd###.mat
	 % store as 32 bits sound (RpvsdEx SerSource uses 32 bit words,
	 % although RP2 has 24 bit AD convertor) 
% 	audiowrite(filename,snd,round(Fs),'BitsPerSample',24);
	soundFs = round(Fs);
	save(filename,'snd','soundFs');
	disp(filename)
	
	plot(snd)
	ylim([-12 12])
	drawnow
	end
end

% [snd,Fs] = audioread(filename);
% isfloat(snd)
% plot(snd)
% ylim([-1 1])
p		= audioplayer(snd,soundFs);
playblocking(p);

