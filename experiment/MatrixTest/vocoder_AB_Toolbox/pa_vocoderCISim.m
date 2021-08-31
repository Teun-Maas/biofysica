function [y,slowEnv,slowEnvE] = pa_vocoderCISim(x,varargin)
% [y,slowEnv,slowEnvE] = vocoderCISim(x, {'params',val,...})
%
% Inputs
%   x           Input sound; assumed sampling rate of 17400 Hz.
% Parameters
% === Analysis
%   nChans      Number of channels (default: 15)
%
% === Reconstruction
%   reconStyle  Style to use for reconstruction filter. 0: trapezoidal, 1: triangular, 2: exponential.
%               (default: 1)
%   reconSlope  Slope (dB/oct) (default: 40)
%   bSyncPhase  If true, synchronize the phase to the sinusoidal wave
%   (default: 0)
%   bSynthSin   Synthesize with sinusoid (default: 0).
%   nSteps      Number of steps to digitize to in each channel.
%   IDR         Dynamic range of the strategy. Determines number of steps.
%   bRemoveBelowIDR  Set values below IDR to very low value. Otherwise, set
%   at -IDR.
% Addition of effect of limited dynamic range and number of steps. This
% would allow us to study the effect of that variable on performance.
% Scaling:
% Ideally, want to do like the processor. We don't want to include  AGC.
% - For full-scale input: 60 dB SPL (AGC knee).
%   - Put that at the top of the IDR range.
%   - Instantaneously, therefore, we will be able to go above IDR top.
%   Limit that by 12 dB. Assume that the steps in that range are
%   logarithmic.
% - Inputs greater than full-scale RMS are reduced to full scale.

% (c) 2007 Litvak, Advanced Bionics
% 2012 Modified by: Marc van Wanrooij

%% Analysis parameters
%% Optional arguments
nChans       = keyval('nChans',varargin); % Number of CI channels
if isempty(nChans)
	nChans		= 15; % Default = 15, MW: in practice?
end

%% Reconstruction parameters
reconStyle	= 1; % triangular
reconSlope	= 40; % Decay (dB/oct)
bSyncPhase	= 0;
bSynthSin	= 0;

%% Compression
nSteps				= 100;
IDR					= 80; % Input dynamic range
sigLevelAtIDRTop	= 10^(-27.3/20); % Level of signal (rms) at top of IDR range.
powerCoef			= 1;
bRemoveBelowIDR		= 1; % Flag to remove power below IDR

%% Parse the parameters
vArgIn = varargin;
for ind = 1:(length(vArgIn )/2)
	if(~exist(vArgIn {ind*2-1},'var'))
		error(['Optional parameter ' vArgIn{ind*2-1} ' doesn''t exist']);
	end
	eval([vArgIn{ind*2-1} '= varargin{ind*2} ;']);
end

%% Parameters which are not settable
nFFT		= 256;   % Size of the FFT
hf			= 0.25;
fLow		= 200; % (Hz)
fHigh		= 10000; % (Hz)
bLimitLast	= 0;
Fs			= 17400;

% figure(668)
% subplot(211)
% specplot(x,Fs);

%% Determine division into channels (if necessary)
f						= spacing_log(fLow,fHigh,nChans);
fL						= f(1:nChans);
fH						= f(2:nChans+1);
[~,indxStart,indxEnd]	= pa_computeClosestFFTFiltGl(fL, fH, Fs/nFFT);
% This fix will help with small number of channels
% iStartArray(1) = round(fLow / (FS/nFFT));

% Extend the end array (by default to 8000)
symEndArray				= max(indxEnd);
mxInd					= nFFT;
indxEnd(nChans)			= nFFT/2-3;

sFitParams.nChans		= nChans;
sFitParams.fLow			= fLow;
sFitParams.fHigh		= fHigh;
sFitParams.bLimitLast	= bLimitLast;
sFitParams.FS			= Fs;
sFitParams.nFFT			= nFFT;

sFitParams.iStartArray	= indxStart;
sFitParams.iEndArray	= indxEnd;
sFitParams.mxInd		= mxInd;
sFitParams.symEndArray	= symEndArray;
bFitUpdated				= 1;

indxStart				= sFitParams.iStartArray;
indxEnd					= sFitParams.iEndArray;
mxInd					= sFitParams.mxInd;

% With hann windowing on both input and output,
% we need 25% window overlap for smooth reconstruction
hop						= nFFT*hf;
% Effect of hanns at both ends is a cumulated cos^2 window (for
% r = 1 anyway); need to scale magnitudes by 2/3 for
% identity input/output
scf						= 2/3;

% Calculate the basic STFT, magnitude scaled
X						= scf*stft(x(:)',nFFT,nFFT,hop);
nT = size(X,2);
%X = X.^2;

%% Divide into channels
%mxInd = 1000;
slowEnv = NaN(nChans,nT);
for iChan = 1:nChans
	iEnd = min([mxInd indxEnd(iChan)]);  % TODO: fix that should be put into code
	if(indxStart(iChan) == indxEnd(iChan))
		slowEnv(iChan,:) = 0.5*log2((abs(X(indxStart(iChan),:)).^2));
	else
		slowEnv(iChan,:) = 0.5*log2(sum(abs(X(indxStart(iChan):iEnd,:)).^2));
	end
end

% Below is full scale empirically determined. This is done by passing a
% tone through simulation at full scale, and looking and slowEnv level.
IDRtop		= 5.705 + log2(sigLevelAtIDRTop); %log2(nFFT*2/3/4);
IDRlog2		= (IDR/20)*log2(10);
stepSize	= IDRlog2/(nSteps-1);

%% Digitize each to the nearest step.
iBelow			= find(slowEnv < IDRtop-IDRlog2);
slowEnv(iBelow) = IDRtop-IDRlog2;
slowEnv(slowEnv > IDRtop+2) = IDRtop+2;

%% Compress slow envelope
slowEnvNorm		= (slowEnv - (IDRtop-IDRlog2)) / (IDRlog2);
slowEnvNorm		= slowEnvNorm.^powerCoef;
slowEnv			= slowEnvNorm*IDRlog2 + (IDRtop-IDRlog2);

%% Determine maximum step given quantization, and fix slowEnv to that.
% with bug in still: slowEnv			= interp1(IDRtop-IDRlog2:stepSize:IDRtop+2,IDRtop-IDRlog2:stepSize:IDRtop+2, slowEnv,'nearest');
slowEnv         = interp1(IDRtop-IDRlog2:stepSize:IDRtop+2+stepSize,IDRtop-IDRlog2:stepSize:IDRtop+2+stepSize, slowEnv,'nearest');
%% Remove steps below IDR.
if(bRemoveBelowIDR)
	slowEnv(iBelow) = -inf;
end

%% Replace 0s by -80
for iChan = 1:nChans
	%     iInf = find(isinf(slowEnv(iChan,:)));
	iInf = find(isinf(slowEnv(iChan,:)) | isnan(slowEnv(iChan,:))); %MW
	
	slowEnv(iChan,iInf) = -80 * ones(size(iInf));
end

slowEnvE = slowEnv;

%% Now resynthesize the waveform

sReconParams.reconStyle = -1;
sReconParams.reconSlope = -1;


%% Create reconstruction filters
if(bFitUpdated || sReconParams.reconStyle ~= reconStyle || ...
		sReconParams.reconSlope ~= reconSlope )
	% Make sure that last expansion is symmetric
	iEnd = indxEnd;
	iEnd(length(iEnd)) = sFitParams.symEndArray;
	sReconParams.reconFilters = createReconFilters(reconStyle,reconSlope,nFFT,Fs,indxStart,iEnd);
	sReconParams.reconStyle = reconStyle;
	sReconParams.reconSlope = reconSlope;
end

%% Reconstruct the envelope
Y = zeros(nFFT/2+1, size(slowEnv,2));
env = 2.^(slowEnvE);

for iChan = 1:nChans
	Y = Y + sReconParams.reconFilters(:,iChan)*env(iChan,:);
end

if bSynthSin
	% Average the energy in the bands; store it in the middle of the bands
	Ysin = zeros(size(Y));
	
	for iChan = 1:nChans
		iEnd = min([mxInd indxEnd(iChan)]);  % TODO: fix that should be put into code
		iMid = floor((indxStart(iChan)+iEnd)/2);
		
		if(indxStart(iChan) == iEnd)
			Ysin(iMid,:) = sqrt(abs(Y(indxStart(iChan),:)).^2);
		else
			Ysin(iMid,:) = sqrt(sum(abs(Y(indxStart(iChan):iEnd,:)).^2));
		end
	end
	
	Y = Ysin;
end
if bSyncPhase
	t = ((1:size(Y,2))-1)*nFFT/4/Fs;
	freqs = -(0:nFFT/2) * Fs/nFFT;
	phase = 2*pi*freqs(:)*t;
	
else
	% Add random phase
	phase = rand(size(Y))*2*pi;
end
Y = Y .* exp(1i*phase);


%% Invert to a waveform
y = istft(Y, nFFT, nFFT, hop)';

%% PLOT THIS
% figure(668)
% subplot(212)
% specplot(y,Fs)

function specplot(y,Fs)

% specgram(y,256,Fs);

snd			= y;
nsamples	= length(snd);
t			= nsamples/Fs*1000;
dt			= 12.5;
nseg		= t/dt;
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 1000;
% figure;
spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
% colorbar
cax = caxis;
caxis([0.7*cax(1) 1.1*cax(2)])
ylim([100 8000])
set(gca,'YTick',[0.5 1 2 4 8 16]*1000);
set(gca,'YTickLabel',[0.5 1 2 4 8 16]);
set(gca,'YScale','log');