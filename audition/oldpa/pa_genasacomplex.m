function [S,M,Fs] = pa_genasacomplex(stim,Lvl)
% [S,M,FS] = PA_GENASACOMPLEX(STIM,LVL)
%
if nargin<1
	stim = 4;
end
if nargin<2
	Lvl = 1;
end
%% Creating stimulus containing signal and masker
% This function generates an synchronous informational Multiple Bursts Same
% for a given level (Lvl).
% S is the sound containing signal and maskers
% M is the sound containing only maskers

%%  Setting stimulus sound parameters
dur			= 20; % ms
Fs          = 48828.125; % TDT Nyquist sampling frequency (Hz)
N			= round(dur/1000*Fs);  % samples
durenv		= 2.5; %ms
Nenv		= round(durenv/1000*Fs); %samples envelope

%% Creating frequency bands, 0.2-14 kHz spaced 1/6 octave
F1			= 200;					%beginning frequency hz
Fsig		= 4000;					%frequency signal (Hz)
Nfreq		= 38;					%number of frequency bands
Nmask		= 4;					%number of tones in masker
freqband	= pa_oct2bw(F1,(0:Nfreq-1)*1/6);	%frequency band, 0.2-14 kHz spaced 1/6 octave
protband	= pa_oct2bw(Fsig,[-1/3 1/3]);		%create protected band (1/3 octave above/below Fsig)
sel			= freqband<protband(1) | freqband>protband(2);
freqband	= freqband(sel);
Nfreq		= numel(freqband);			%Nfreq - protected band

%% Creating temporal bands
rate		= 10;								% cycli/s frequency
period		= 1/rate;							% s time interval
Nrep		= 4;								% repetitions of signal
Nperiod		= round(period*Fs);					% time interval between signals (samples)

% timebin		= floor(nperiod/N); % samples

%% Timings and frequencies
onset	= 0;
T		= pa_rndval(1,Nperiod,[Nmask,Nrep]);
F		= T;
for ii	= 1:Nrep
	indx	= randperm(Nfreq,Nmask);	%selecting (4) maskers from Nfreq
	F(:,ii) =  sort(freqband(indx));		%retrieving corresponding frequencies from indx (masker frequencies) [Hz]
end

%% Synchronous
switch stim
	case 1
	case 2
		F(:,2:Nrep) = repmat(F(:,1),1,(Nrep-1));
		T(:,2:Nrep) = repmat(T(:,1),1,(Nrep-1));
	case 3
		T = zeros(Nmask,Nrep);
	case 4
		F(:,2:Nrep) = repmat(F(:,1),1,(Nrep-1));
		T			= zeros(Nmask,Nrep);
end

%% Repetitive Signal
S		= zeros(1,round((Nrep+1)*Nperiod));			%creating empty vector for sounds
M		= zeros(1,round((Nrep+1)*Nperiod));			%creating empty vector for maskers
s		= pa_gentone(Fsig,dur,'nenvelope',Nenv,'Fs',Fs);		%#ok<*NASGU> %signal synchronous iMBS
for ii = 1:Nrep
	indx = ((ii-1)*Nperiod+1):((ii-1)*Nperiod+N);
	S(indx) = Lvl*s;
end

%% Masker
for ii = 1:Nrep
	for jj = 1:Nmask
		s = pa_gentone(F(jj,ii),dur,'nenvelope',Nenv,'Fs',Fs);		%#ok<*NASGU> %signal synchronous iMBS
		indx = ((ii-1)*Nperiod+1+T(jj,ii)):((ii-1)*Nperiod+N+T(jj,ii));
		S(indx) = S(indx)+s;
		M(indx) = M(indx)+s;
		
	end
end


pa_plotspec(S,Fs)
hold on
pa_horline(protband,'w-');
drawnow;

wavplay(M,Fs);
wavplay(S,Fs);

