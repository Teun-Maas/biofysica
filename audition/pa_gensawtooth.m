function [x,Fs] = pa_gensawtooth(k,f0)

Fs	= 44100;

t = 0:1/Fs:.5;
if nargin<1
	k = 5;
end
if nargin<2
f0 = 100;
end

for ii = k
	figure(1)
	subplot(211)
	cla
	x = sawtooth(ii,f0,t);
	
	N = round(100/1000*Fs);
	x = pa_envelope(x',N);
	
	plot(t,x);
	xlim([0 0.5]);
	ylim([-1.3 1.3]);

	P = specgram(x,Fs);
	drawnow
end
% f0cor
	p = audioplayer(x,Fs);
	playblocking(p);

% 	y = pa_psola(x
function x = sawtooth(k,f,t)

x = 0;
for ii = 1:k
	x = x -2/pi * sin(2*ii*f*t)/ii;
end


function P = specgram(Y,Fs)
di = 1;
dt = 1;
nsamples	= length(Y);
t			= nsamples/Fs*1000;
nseg		= round(t*di/dt)
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 2^13;

[S,F,T] = spectrogram(Y,window,noverlap,nfft, Fs);
P	= S.*conj(S);
P = log10(P+1);

subplot(212)
cla
	mxscal = max(abs(P(:)));

	imagesc(T,F,P);
	xlabel('Time (s)');
	ylabel('Frequency (Hz)');
	caxis([-mxscal mxscal]);
	xlim([0 0.5]);
	ylim([50 3000]);
	set(gca,'YDir','normal');
	