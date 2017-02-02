function [Mag,Freq,Azimuth,Elevation] = sweep2mag(Sweep,Wav,varargin)

% [MAG,FREQ,AZIMUTH,ELEVATION] = SWEEP2MAG(SWEEP,WAV)
%
% See also READHRTF, GENSWEEP


%% Initialization
chan		= keyval('channel',varargin,1); % channel
NFFT		= keyval('nfft',varargin,1024); % Number of FFT samples/bins (each sweep should contain NFFT samples)
NSweep		= keyval('nsweep',varargin,18); % the number of presented sweeps -2 (beginning and end are removed, as these contain on and offset ramps)
Fs			= keyval('Fs',varargin,48828.125); % sampling frequency (of sound and of recording)
nBegin		= keyval('nbegin',varargin,1); % begin sample
GraphFlag	= keyval('display',varargin,false); % visualize
Azimuth     = keyval('azimuth',varargin);
Elevation   = keyval('elevation',varargin);
maxLag		= keyval('maxLag',varargin,1000); % maximum lag
Freq		= (0:(NFFT-1))*Fs/NFFT;
nloc		= size(Sweep,3); % number of presentations, typically 1 per location
nchan		= size(Sweep,2);
thresh		= 4;

%% Alignment 
% First, we check by crosscorrelation whether there is a delay of the sound
% present in the recording. We need this, as we want to analyze only that
% part of the recording that contains sweeps.

L       = NaN(1,nloc);

% for each location locIdx
for locIdx = 1:nloc
	if nchan>1
		d           = squeeze(Sweep(:,chan,locIdx)); % obtain the measured Sweep data from channel chan
	else
		d           = squeeze(Sweep(:,locIdx)); % obtain the measured Sweep data from channel chan
	end
	[c,lags]    = xcorr(d,Wav,maxLag,'none');
	[~,indx]    = max(abs(abs(c))); % maximum crosscorrelation
	lag         = lags(indx); % determines lag @ max(cc)
	L(locIdx)	= lag;
	
% 	graph1(d,Wav,lag,lags,c,Fs,locIdx,nBegin,NFFT,GraphFlag,600)
end
L			= L(L<500);
lag			= round(nanmedian(L));

%% FFT & Reshaping
% We now reshape the recording data, so that we can use matrix
% manipulation. We will create NFFT by Nsweep large matrices.
% Determine Magnitude from the time-domain data
% For each sweep, determine spectrum via fft
% THEN average over sweeps
Mag			= NaN(NFFT,nloc);
for locIdx = 1:nloc
	if nchan>1
		d           = squeeze(Sweep(:,chan,locIdx)); % obtain the measured Sweep data from channel chan
	else
		d           = squeeze(Sweep(:,locIdx)); % obtain the measured Sweep data from channel chan
	end
	nOnset      = lag+nBegin; % Onset of sound in the recording in samples
	indx        = nOnset + NFFT + (1:(NFFT*NSweep)); % sample index: start at sound, +1 sweep, until Nsweep
	data        = d(indx); 
	data        = reshape(data, NFFT, NSweep); % reshaping the data so that it has (usually) 1024 samples and (about) 18 sweeps
	
	graph1(data(:),Wav,0,lags,c,Fs,locIdx,nBegin,NFFT,GraphFlag,601)

%% Correct for errors in Time
% 	% We might need to correct for noises (e.g. subject swallowing)
% 	meandata	= mean(data,2); % mean over sweeps
% 	stddata		= std(data,0,2); % std over sweeps
% 	deltadata	= bsxfun(@minus,data,meandata);
% 	
% 	sel			= (abs(deltadata) - thresh*repmat(stddata,1,size(data,2))) >0; % remove sweep when any part exceeds the threshold
% 	sel			= ~logical(sum(sel));
% 	data		= data(:,sel);


	%% Spectrum
	Spec		= NaN(size(data));
	x			= bsxfun(@minus,data',mean(data)')'; % removes DC, offset, F=0 component
	nfft        = 2^(nextpow2(length(x))); % should be 1024
	for jj		= 1:size(data,2)
		Spec(:,jj)  = abs(fft(x(:,jj),nfft));
	end
	Mag(:,locIdx)   = mean(Spec,2);
	

end
% NumUniquePts	= ceil((NFFT+1)/2);
% Spec			= Spec(1:NumUniquePts,:);

%% Average over location repetition
Locations	= [Azimuth Elevation];
uLoc		= unique(Locations,'rows');
M			= NaN(size(Mag,1),size(uLoc,1));
muall		= mean(Mag,2);
for ii = 1:size(uLoc,1)
	sel		= Azimuth == uLoc(ii,1) & Elevation == uLoc(ii,2);
	m		= Mag(:,sel);
	mu		= mean(m,2);
	M(:,ii) = mu;
	
	if GraphFlag==2
		MUall	= repmat(muall,1,size(m,2));
		m	= m./MUall;
		m	= hsmooth(m);
		mu2 = mean(m,2);
		figure(603)
		clf
		subplot(121)
		plot(Freq(1:513),m,'k-','Color',[.7 .7 .7]);
		hold on
		plot(Freq(1:513),mu2,'k-');
		x = [1000 2000 4000 8000 16000];
		set(gca,'Xtick',x,'XTickLabel',x/1000);
		axis square
		xlim(minmax(x));
		
		subplot(122)
		mu2		= repmat(mu2,1,size(m,2));
		m		= m-mu2;
		plot(Freq(1:513),m,'k-','Color',[.7 .7 .7]);
		hold on
		x = [1000 2000 4000 8000 16000];
		set(gca,'Xtick',x,'XTickLabel',x/1000);
		xlim(minmax(x));
		ylim([-4 4])
		horline
		axis square
% 		pause
	end
	
end

%% Output
Mag			= M;
Azimuth		= uLoc(:,1);
Elevation	= uLoc(:,2);

function graph1(d,Wav,lag,lags,c,Fs,locIdx,nBegin,NFFT,GraphFlag,fignr)
if GraphFlag
	% first some scaling
	d = d./max(d(:));
	Wav = Wav./max(Wav(:));
	td = 0:length(d)-1; td = 1000*(td-lag)/Fs;
	tw = 0:length(Wav)-1; tw = 1000*tw/Fs;
	figure(fignr)
	clf
	subplot(211)
	plot(tw,Wav,'r-');
	hold on
	plot(td,d,'k-');
	xlim([0 max(td)]);
	ylim([-3 3]);
	horline(rms(d))
	title(locIdx)
	verline(nBegin+NFFT);
	verline(nBegin+2*NFFT);
	
	subplot(212)
	plot(lags,c);
	drawnow
end