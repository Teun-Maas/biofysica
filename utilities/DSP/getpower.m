function [f,mx,ph,a,h]=getpower(x,Fs,varargin)
% [F,A,PH] = PA_GETPOWER(X,FS)
%
% Get power A and phase PH spectrum of signal X, sampled at FS Hz.
%
% PA_GETPOWER(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'display'	- display graph. Choices are:
%					0	- no graph (default)
%					>0	- graph
%	'color'	- specify colour of graph. Colour choices are the same as for
%	PLOT (default: k - black).


% 2011  Modified from Mathworks Support:
% http://www.mathworks.com/support/tech-notes/1700/1702.html
% by: Marc van Wanrooij

%% Initialization
if nargin<2
	Fs = 48828.125;
end
% Optional display arguments
disp		= keyval('display',varargin,false);
orient		= keyval('orientation',varargin,'x');
col			= keyval('color',varargin,'k');


%%
% Time vector of 1 second
% Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
nfft			= keyval('nfft',varargin,2^(nextpow2(length(x))));
NumUniquePts	= ceil((nfft+1)/2);

maxnfft			= 2^(nextpow2(length(x))-1);
nrep			= floor(maxnfft/nfft);
if nrep>1
	idx = (nfft*nrep);
	x = x(1:idx);
	
	x		= reshape(x,nfft,nrep);
	
	%%
	fftx = 0;
	for ii = 1:nrep
		% 		[mx,ph] = getspectrum(x,nfft,NumUniquePts);
		% Take fft, padding with zeros so that length(fftx) is equal to nfft
		fftx			= sqrt(fftx.^2+fft(x,nfft).^2);
		% Calculate the number of unique points
		% FFT is symmetric, throw away second half
		% Take the magnitude of fft of x and scale the fft so that it is not a function of
		% the length of x
	end
	
	fftx			= fftx(1:NumUniquePts)/nrep;
	
	mx				= abs(fftx)/length(x);
	ph				= angle(fftx);
	
	% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
	% The DC component and Nyquist component, if it exists, are unique and should not
	% be mulitplied by 2.
	if rem(nfft, 2) % odd nfft excludes Nyquist point
		mx(2:end)	= mx(2:end)*2;
	else
		mx(2:end -1) = mx(2:end -1)*2;
	end
	
	%%
	% 	keyboard
else
	[mx,ph] = getspectrum(x,nfft,NumUniquePts);
end
% This is an evenly spaced frequency vector with NumUniquePts points.
f		= (0:NumUniquePts-1)*Fs/nfft;
a = mx;
% Take the square of the magnitude of fft of x -> magnitude 2 power
% mx				= mx.^2;
mx		= 20*log10(mx);
sel		= isinf(mx);
mx(sel) = min(mx(~sel));

%% Display option
if disp
	if strcmpi(orient,'x')
		% 		w = aweight(f);
		% 		m = mx+w;
		m = mx;
		h = semilogx(f,m,'-');
		set(h,'Color',col);
		hold on
		set(gca,'XTick',[0.05 0.5 1 2 3 4 6 8 10 14]*1000,...
			'XTickLabel',[0.05 0.5 1 2 3 4 6 8 10 14]);
		title('Power Spectrum');
		xlabel('Frequency (kHz)');
		ylabel('Power (dB)');
		nicegraph;
		
		%%
	elseif strcmpi(orient,'y')
		h = semilogy(mx,f);
		set(h,'Color',col);
		% 		set(gca,'YTick',[0.05 1 2 3 4 6 8 10 14]*Fs);
		% 		set(gca,'YTickLabel',[0.05 1 2 3 4 6 8 10 14]);
		title('Power Spectrum');
		ylabel('Frequency (Hz)');
		xlabel('Power (dB)');
		% 	axis square;
	end
end

function [mx,ph] = getspectrum(x,nfft,NumUniquePts)
% Take fft, padding with zeros so that length(fftx) is equal to nfft
fftx			= fft(x,nfft);
% Calculate the number of unique points
% FFT is symmetric, throw away second half
fftx			= fftx(1:NumUniquePts);
% Take the magnitude of fft of x and scale the fft so that it is not a function of
% the length of x
mx				= abs(fftx)/length(x);
ph				= angle(fftx);

% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not
% be mulitplied by 2.
if rem(nfft, 2) % odd nfft excludes Nyquist point
	mx(2:end)	= mx(2:end)*2;
else
	mx(2:end -1) = mx(2:end -1)*2;
end