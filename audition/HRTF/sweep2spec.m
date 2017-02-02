function [Spec,Freq,Tijd] = pa_sweep2spec(Sweep,Wav,chan,NFFT,NSweep,Fs,GraphFlag)
% [SPEC,FREQ,TIME] = SWEEP2SPEC(SWEEP,WAV,CHAN,NFFT,NSWEEP,CHAN)
%
% See also READHRTF, GENSWEEP


%% Initialization
if nargin<3
	chan = 1;
end
if nargin<4
	NFFT    = 1024;
end
if nargin<5
	NSweep = 18;
end
if nargin<6
	Fs     = 48828.125;
end
if nargin<7
	GraphFlag = 0;
end

nloc	= size(Sweep,3);
nchan	= size(Sweep,2);

%% FFT & Reshaping
nBegin  = ceil(20/1000*Fs);
Tijd    = NaN(NFFT,nloc);
Spec    = Tijd;
L       = NaN(1,nloc);

% for each location i
for i = 1:nloc
	if nchan>1
		d           = squeeze(Sweep(:,chan,i)); % obtain the measured Sweep data from channel chan
	else
		d           = squeeze(Sweep(:,i)); % obtain the measured Sweep data from channel chan
	end
	d           = d(101:end); % Remove the first 100 samples, somehow HumanV1 adds 100 extra samples
	% Rough Alignment Method 1
	if length(d)==length(Wav)
		[c,lags]    = xcorr(d,Wav,'coeff');
	else
		[c,lags]    = xcorr(d,Wav,'none');
	end
	[~,indx]    = max(abs(c));
	lag         = lags(indx);
	
	if GraphFlag
		td = 0:length(d)-1; td = 1000*(td-lag)/Fs;
		tw = 0:length(Wav)-1; tw = 1000*tw/Fs;
		figure(600)
		clf
		subplot(211)
		plot(tw,Wav,'r-');
		hold on
		plot(td,d,'k-');
		xlim([0 max(td)]);
		ylim([-3 3]);
		subplot(212)
		plot(lags,c,'k-');
		xlim([min(lags) max(lags)]);

		drawnow
		pause;
		
	end
	L(i)        = lag;	
end
if GraphFlag
	L = L(L<250);
	%% Check whether lags are correct
	mnL = min(L);
	mxL = max(L);
	stp = std(L)./sqrt(numel(L));
	x = unique(round(min(L):stp:max(L)));
	figure(601)
	N = hist(L,x);
	[~,indx] = max(N);
	bar(x,N,1);
	verline(round(median(L)));
	verline(x(indx));
end
lag = median(L);
% Aligmnent Method 2
% If Sweep is too weak, the time delay might be judged incorrectly
% Correct for this:
if lag>250 || lag<1
	lag     = 210;
end
for i = 1:nloc
	if nchan>1
		d           = squeeze(Sweep(:,chan,i)); % obtain the measured Sweep data from channel chan
	else
		d           = squeeze(Sweep(:,i)); % obtain the measured Sweep data from channel chan
	end
	d           = d(101:end); % Remove the first 100 samples, somehow HumanV1 adds 100 extra samples
	nOnset      = lag+nBegin;
	indx        = nOnset + NFFT + (1:(NFFT*NSweep));
	data        = d(indx);
	data        = reshape(data, NFFT, NSweep);
	meandata    = mean(data,2)';
	meandata    = meandata - mean(meandata);
	Tijd(:,i)   = meandata;
	nfft        = 2^(nextpow2(length(meandata)));
	X           = fft(meandata,nfft);
	Spec(:,i)   = X;
end
% NumUniquePts	= ceil((NFFT+1)/2);
% Spec			= Spec(1:NumUniquePts,:);
Freq			= (0:(NFFT-1))*Fs/NFFT;

