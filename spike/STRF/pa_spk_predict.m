function [S,C,P] = pa_spk_predict(strf,Y,dt,varargin)
% [S,C,P] = PA_SPK_PREDICT(STRF,Y)
%
% Prediction S of firing rate to stimulus Y of cell with receptive field
% STRF.
%
% Stimulus Y has spectrogram P
% Spectrogram of prediction is C (convolution of P and STRF)
%
% See also PA_SPK_RIPPLE2STRF

% 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
dspFlag       = pa_keyval('display',varargin);  % sample frequency of ripple stimulus
if isempty(dspFlag)
	dspFlag = 1;
end
meth       = pa_keyval('method',varargin);  % sample frequency of ripple stimulus
if isempty(meth)
	meth = 1;
end
hfshift       = pa_keyval('shift',varargin);  % sample frequency of ripple stimulus
if isempty(hfshift)
	hfshift		= 0;
end
Fs       = pa_keyval('Fs',varargin);  % sample frequency of ripple stimulus
if isempty(Fs)
	Fs		= 50000;
end
absFlag       = pa_keyval('abs',varargin);  % sample frequency of ripple stimulus
if isempty(absFlag)
	absFlag		= 0;
end
bf      = pa_keyval('bf',varargin);  % sample frequency of ripple stimulus
if isempty(bf)
	bf		= 1;
end
smth      = pa_keyval('smooth',varargin);  % sample frequency of ripple stimulus
if isempty(smth)
	smth		= 0;
end
di       = pa_keyval('interp',varargin);  % sample frequency of ripple stimulus
if isempty(di)
	di		= 1;
end
nsamples	= length(Y);
t			= nsamples/Fs*1000;
nseg		= t*di/dt;
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 1000;
xoct		= linspace(0.25,2.5,10);

% Old
% window = round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
% noverlap	= 2^0; % temporal smoothing
% noverlap	= 2^(nextpow2(0.6*window)-1)
% nfft		= 2^10;

%% Remove late strf
% strf = strf(:,1:6);
% mn = -1;
% sel = strf<mn;
% strf(sel) = mn;

%% STRF
[strf, kshift]	= pa_spk_strfshiftmax(strf); 
% kshift = 0;
% kshift = kshift-3; % Improve?
f0				= 250*2^(-kshift*0.25+2.5*(bf-1)) ;
f				= pa_oct2bw(f0,xoct);

if di>1
	n			= size(strf,2);
	x			= (1:n)*t;
	y			= xoct;
	z			= strf;
	XI			= linspace(min(x),max(x),n*di);
	YI			= linspace(min(y),max(y),n*di);
	ZI			= interp2(x,y',z,XI,YI','spline');
	
	% rescale
	xoct		= YI;
	strf		= ZI;
	dt			= dt/di;
end

if dspFlag
	t		= (1:size(strf,2))-1;
	t		= t*dt;
	mxscal	= max(abs(strf(:)));
	ytck	= linspace(min(xoct),max(xoct),6);
	ylbl	= round(pa_oct2bw(f0,ytck)/1000*10)/10;
	xtck	= round(linspace(min(t),max(t),6));

	figure(dspFlag)
	subplot(321)
	imagesc(t,xoct,strf);
	set(gca,'YDir','normal','YTick',ytck,'YTickLabel',ylbl,'XTick',xtck);
	axis square;
	ylabel('Frequency (kHz)');
	caxis([-mxscal mxscal]);
	
	subplot(323)
	imagesc(t,xoct,strf);
	set(gca,'YDir','normal','YTick',ytck,'YTickLabel',ylbl,'XTick',xtck);
	axis square;
	ylabel('Frequency (kHz)');
	caxis([-mxscal mxscal]);
end

%% Determine spectrogram of stimulus
% For convolution with STRF we require logarithmically-spaced (octaves)
% frequencies.
[~,~,Ti] = spectrogram(Y,window,noverlap,nfft, Fs);

switch meth
	case 1
		%% Method 1: Interpolate
		[S,F,T] = spectrogram(Y,window,noverlap,nfft, Fs);
		P	= S.*conj(S);
		Pi	= NaN(numel(f),numel(T));
		for ii = 1:size(P,2)
			p = P(:,ii);
			if smth
				p = smooth(p,smth);
			end
			p = interp1(F,p,f,'cubic');
			Pi(:,ii) = p;
		end
		P		= Pi;
	case 2
		%% Method 2: Average NO NO
		[S,F,T]	= spectrogram(Y,window,noverlap,nfft, Fs); %moet geplot worden
		P		= S.*conj(S);
		Pi		= NaN(numel(f),numel(T));
		xoct	= 0.125:0.25:(2.5+0.125);
		fedge	= pa_oct2bw(250,xoct+hfshift);
		for ii = 1:length(fedge)-1
			sel			= F>=fedge(ii) & F<fedge(ii+1);
			p			= P(sel,:);
			p			= nanmean(p);
			Pi(ii,:)	= p;
		end
		P		= Pi;
	case 3
		%% Method 3: Average NO NO
		df = 0.25;
		fWin = [f(1) f(end)];
		[P, T, F] = ic_specgram(Y, fWin, df, dt);
		Pi			= NaN(numel(F),numel(T));
		xoct			= 0.125:0.25:(2.5+0.125);
		fedge		= pa_oct2bw(250,xoct+hfshift);
		for ii = 1:length(fedge)-1
			sel			= F>=fedge(ii) & F<fedge(ii+1);
			p			= P(sel,:);
			p			= nanmean(p);
			Pi(ii,:)	= p;
		end
		P		= Pi;
end
P		= log10(P+1);
% P		= P-min(P(:)); % QUE?

if dspFlag % Plot spectrogram
	mxscal = max((P(:)));
	ytck	= linspace(min(xoct),max(xoct),6);
	ylbl	= round(pa_oct2bw(f0,ytck)/1000*10)/10;

	subplot(322)
	imagesc(T*1000,xoct,P);
	xlabel('Time (ms)'); ylabel('kHz');
	set(gca,'YDir','normal','YTick',ytck,'YTickLabel',ylbl);
	ylabel('Frequency (kHz)');
	colorbar;
	caxis([-mxscal mxscal]);
end

%% Convolution
[m,n]	= size(P);
C		= NaN(m,size(strf,2)+n-1);
for ii = 1:m
	ci		= conv(strf(ii,:),P(ii,:));
	C(ii,:) = ci;
end
S = sum(C); % Prediction
S = S(1:n); % Remove trailing conv artefacts
C = C(:,1:n);

if absFlag
	sel		= S<0;
	S(sel)	= 0;
end

S = S*dt;

if dspFlag
	t		= 0:(length(C)-1);
	t		= t*dt;
	mxscal	= max(abs(C(:)));
	ytck	= linspace(min(xoct),max(xoct),6);
	ylbl	= round(pa_oct2bw(f0,ytck)/1000*10)/10;
	
	subplot(324)
	imagesc(t,xoct,C);
	set(gca,'YDir','normal','YTick',ytck,'YTickLabel',ylbl);
	xlabel('Time (ms)'); ylabel('kHz');
	caxis([-mxscal mxscal]);
	colorbar;
	
	subplot(326)
	plot(t,S,'k-','LineWidth',2);
	pa_horline;
	xlim([min(t) max(t)]);
	colorbar;
end

