function snd = gensweep(Nsample, Fstart, Fstop, Fs, Nsweep, varargin)
% Generate Schroeder Sweep
%
% snd = PA_GETSWEEP(N,FSTART,FSTOP)
%
%   return a sweep of N samples with a flat magnitude
%   spectrum between frequencies FSTART and FSTOP
%
% See also HOW2CREATESOUND, PA_GENGWN, PA_ENVELOPE, PA_WRITEWAV

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<1
	Nsample  = 2^10; %1024 samples
end
if nargin<2
	Fstart  = 0; % Hz
end
if nargin<4
	Fs     = 48828.125;
end
if nargin<3
	Fstop   = Fs; % Hz
end
if nargin<5
	Nsweep = 21;
end
%% Optional arguments
dspFlag       = keyval('display',varargin);
if isempty(dspFlag)
	dspFlag	= 0;
end
plee       = keyval('play',varargin);
if isempty(plee)
	plee	= 'n';
end

%% Obtain single sweep
% patchwork: Lowpeak generates a sweep of 2N samples for
%            a magitude spectrum of N samples
N           = round(0.5*Nsample);
% OK, now let's fill that magnitude spectrum
Magnitude   = zeros(1,N);
nstart      = max(1,round(N*Fstart/Fs));
nstop       = min(N,round(N*Fstop/Fs));
Magnitude(nstart:nstop) = ones(1,(nstop-nstart+1));
snd           = lowpeak(Magnitude);

%% Process Sweep
snd = repmat(snd,1,Nsweep);

%% Graphics
if dspFlag
	close all
	figure
	t = (1:length(snd))/Fs*1000;
	subplot(221)
	plot(t,snd,'k-')
	ylabel('Amplitude (au)');
	ylabel('Time (ms)');
	xlim([min(t) max(t)]);
	
	subplot(224)
	pa_getpower(snd,Fs,'orientation','y');
	ylim([0 20000])
	set(gca,'Yscale','linear','YTick',(0:5:20)*1000,'YTickLabel',0:5:20);
	ylabel('Frequency (kHz)');
	
	subplot(223)
	nfft		= 2^11;
	window		= 2^6; % resolution
	noverlap	= 2^4; % smoothing
	spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
	set(gca,'Yscale','linear','YTick',(0:5:20)*1000,'YTickLabel',0:5:20);
	
	ylim([0 20000])
	ylabel('Frequency (kHz)');
	
	% 	set(gca,'YTick',[0.05 1 2 3 4 6 8 10 14]*1000);
	% 	set(gca,'YTickLabel',[0.05 1 2 3 4 6 8 10 14]);
	
	drawnow
	
end

%% Play
if strcmpi(plee,'y');
	sndplay = pa_envelope(snd',round(10*Fs/1000));
	p		= audioplayer(sndplay,Fs);
	playblocking(p);
end
function Signal = lowpeak(Magnitude)
%  FUNCTION X = LOWPEAK (M)
%
%  DESCRIPTION
%    Returns the time signal X with minimized peakfactor
%    given a frequency magnitude response M. The algorithm is
%    extracted from Schroeder et al (1970).
%
%  ARGUMENT
%    M  - frequency magnitude response containing the frequency
%         response in equidistant bins which correspond to frequency
%         0 to the Nyquist frequency (best to make the length of M
%         a power of 2). The algorithm calculates a frequency phase
%         response that together with M, yields a low-peak time signal.
%
%  RETURN VALUE
%    X  - Time signal with minimized peak factor
%
%  EXAMPLE
%    >> sweep = lowpeak (ones(1,512));
%
%  .. Hofware Inc...
%
Nbin       = length(Magnitude);
TotalPower = sum(Magnitude.^2);
NormFactor = 1.0/TotalPower;
TwoPi      = 2*pi;
Phi        = 0.0;
Power      = 0.0;
Spectrum   = zeros (1, Nbin);
for j=1:Nbin
	Spectrum(j) = Magnitude(j) * exp (1i*Phi);
	Power = Power + NormFactor*Magnitude(j).^2;
	Phi   = Phi - TwoPi*Power;
end;
Spectrum = [Spectrum -conj([Spectrum(1) Spectrum((Nbin):-1:2)])];
Signal   = imag(ifft(Spectrum));
