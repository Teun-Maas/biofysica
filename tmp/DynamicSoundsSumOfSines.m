function DynamicSoundsSumOfSines
close all
clearvars;

f0 = 0.05;
h = [2 3 7 13 21];

freq	= f0*h;
nfreq	= numel(freq);

[F1,F2] = meshgrid(freq,freq);


xf = unique(F1(:)-F2(:));
xf = xf(xf>0);

%%
Fs			= 48828.125; % Hz
dur			= 20;
nsamples	= floor(dur*Fs);
dur			= nsamples/Fs;
t			= linspace(0,dur,nsamples);

for jj = 1:36
	SoS			= 0;
	for ii = 1:nfreq
		phi = rand(1);
		
		A = sin(2*pi*freq(ii)*t+2*pi*phi);
		SoS = SoS+A;
	end
% 	SoS = SoS-mean(SoS);
	SoS = 2*(SoS-min(SoS))/(max(SoS)-min(SoS))-1;
% 	SoS = SoS/2;
	SoS = ramp(SoS,floor(2*Fs));
% 	wavwrite(SoS,Fs,['snd' num2str(jj) '.wav']);
	
	
	figure(1)
	subplot(211)
	plot(t,SoS);
	hold on

	%%
	v = [0; diff(SoS*45)]*Fs;
	subplot(212)
	plot(t,v);
	hold on

	
% 	pause
	%%
end

figure(1)
nicegraph


%%



%%
% durz = 2;
% nzsamples = floor(durz*Fs);
% z = zeros(1,nzsamples);
% 
% dur = (nsamples+nzsamples)/Fs;
% t = linspace(1,dur,nsamples+nzsamples);

% whos z SoS
% SoS = [z SoS];


%% Generate sound
%Generate GWN with a time ramp function
GWN=gengwn(20);
GWN=ramp(GWN);

%Guardamos el sonido
% wavwrite(GWN,sr,'snd440.wav');
writewav(GWN,'snd440.wav');



function writewav(Signal,fname,ScalingFactor,Fsample,Nbits)
% Write a stimulus to WAV-file
%
% WRITEWAV(SIGNAL,WAVFILE,SF)
%
% Write a stimulus to WAV-file for use with FART
% SF = ScalingFactor
%
% See also GENGWN, GENSWEEP
%

% Copyright 2007
% Author MarcW

%% Initialization
if nargin<3
    ScalingFactor = 0.99;
end
if nargin<4
%   Fsample = 24414.065;
  Fsample = 48828.125;  % TDT Nyquist sampling frequency (Hz)
end;
if nargin<5
    Nbits = 16; % TDT maximal bit depth for WAVs
end
fname = fcheckext(fname,'wav');

%% Minmax signal for amplifier
if ~isempty(ScalingFactor)
    Signal = ScalingFactor.*Signal./max(max(abs(Signal)));
end

%% Write
wavwrite(Signal,Fsample,Nbits,fname);