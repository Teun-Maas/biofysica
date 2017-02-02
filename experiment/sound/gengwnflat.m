function [snd,Fs] = gengwnflat(N, order, Fc, Fn, Fh, varargin)
% function [snd,Fs] = gengwnflat(N, order, Fc, Fn, Fh, varargin)
% Generate Gaussian White Noise Stimulus by defining a flat Magnitude
% spectrum and a random phase.
%
%
% See also GENGWN, WRITEWAV
%

% 2007 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com


%% Initialization

if nargin<5
	Fh          = 500;
end
if nargin<4
	Fn          = 48828.125/2; % TDT Nyquist sampling frequency (Hz)
end
if nargin<3
	Fc          = 20000; % Hz
	Fc			= min([Fc Fn]);
end
if nargin<2
	order       = 500; % samples
end
if nargin<1
	N           = 0.15*Fn*2; % samples
end
Fs = Fn*2;

%% Optional arguments
dspFlag		= keyval('display',varargin,false);
plee		= keyval('play',varargin,'n');

%% Create and Modify Signal
NFFT	= 2^(nextpow2(N));
M		= repmat(100,NFFT/2,1);
M		= [M;flipud(M)];
P		= (rand(NFFT/2,1)-0.5)*2*pi;
P		= [P;flipud(P)];

R		= M.*cos(P);
I		= M.*sin(P);
S		= complex(R,I);

snd		= ifft(S,'symmetric');

% Low-pass filter
snd             = lowpass(snd,'Fc',Fc,'Fs',Fs,'order', order);
% High-pass filter
snd             = highpass(snd,'Fc',Fh,'Fs',Fs,'order', order);

%% Optional Graphics
if dspFlag
	figure;
	disp('>> GENGWN <<');
	subplot(211);
	plot(snd);
	xlabel('Sample number');
	ylabel('Amplitude (a.u.)');

	
	subplot(212)
	getpower(snd,Fn*2,'display',1);
	xlim([50 22000])
end

%% Play
if strcmpi(plee,'y')
	sndplay = envelope(snd',round(10*Fs/1000));
	p		= audioplayer(sndplay,Fs);
	playblocking(p);
end
