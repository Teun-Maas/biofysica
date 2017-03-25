function [newspec, tnew, f] = ic_specgram(s_t, fWin, df, dt)
%	[newspec, tnew, f] = ic_specgram(s_t, f_win, f_bin, t_bin)
%

% Huib Versnel/John van Opstal/Marcel Zwiers
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<4
	dt		= 12.5; % (ms) = 1/(velocity step)
end
if nargin<3
	df		= 0.25; % (oct) = 1/(density step)
end
if (nargin < 2)
	fWin		= [250*2^0.25 250*2^2.5]; % 297 - 1414 Hz
end
corfac	= 1; % Correction?
if dt == 6.25
	s_t		= resample(s_t,40000,50000);
	dt		= 5;
	corfac	= 1.25;
end

halfbin = dt/2;
nshift	= round(50*halfbin); % assume 50 kHz sampling rate
s_new	= [zeros(nshift,1); s_t];

[newspec,tnew,f]	= specgram2(s_new,fWin,df,dt);
tnew				= corfac*tnew;

%% Graphics
if ~nargout
	[logspec,t,f]	= specgram2(s_t, fWin, df, dt);
	t				= corfac*t;
	Nf				= size(logspec,1);
	f1				= (0:Nf-1)*df;
	fticks			= (0:(Nf-1)/4)*df*4;
	flabels			= fWin(1)*2.^fticks;
	
	subplot(2,1,1);
	imagesc(t,f1,log10(logspec));
	set(gca,'YDir','normal','TickDir','out');
	set(gca,'YTick',fticks,'YTickLabel',flabels);
	
	subplot(2,1,2);
	imagesc(tnew,f1,log10(newspec));
	set(gca,'YDir','normal','TickDir','out');
	set(gca,'YTick',fticks,'YTickLabel',flabels);
	set(gcf,'position',[525 40 600 750]);
end

function [logspec, t, f] = specgram2(s_t, fWin, df, dt)
% FUNTION [logspec, t, f] = specgram2(s_t, f_win, f_bin, t_bin)
%
% Average energy in a spectogram with logaritmic frequency bins (overlap in
% time-bins about 1/3; see also specgram). N.B. the averaging over
% frequency-bins is only appropriate for stimuli which are flat in a linear
% frequency-domain.
%
% INPUT:                                                      DEFAULT
%   s_t     - signal in time-domain (sampled at 50 kHz)
%   f_win   - lowest and highest frequency in specgram (Hz)   [250 16000]
%   f_bin   - binsize along frequency-axis (oct)              0.25
%   t_bin   - binsize along time-axis (msec)                  5
%
% OUTPUT:
%   logspec - powermatrix, 't' from left to right, 'f' from bottom to top
%   t       - time
%   f       - frequency
%
% If no output arguments are requested, the dB energyspectrum is plot in a
% figure.

% Marcel Zwiers, 16-1-'01.
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initizalization
if (nargin < 4)
	dt = 5;
end
if (nargin < 3)
	df = 0.25;
end
if (nargin < 2)
	fWin = [250 16000];
end
f1 = fWin(1);
f2 = fWin(2);

% (1 + xover)*t_bin*50 & xover*t_bin*50 => integer (fft)
xover = 0.6;

% calculate the powerspectrogram (about 1/3 overlap)
% [x,nfft,Fs,window,noverlap]
[sig_ft,f,t] = specgram(s_t,(1+xover)*dt*50, 50000,[],xover*dt*50);
pwr_ft       = sig_ft.*conj(sig_ft);

f_cntr		= f1*2.^(0:df:log2(f2/f1));
f_edge		= f1*2.^(-0.5*df : df : log2(f2/f1)+0.5*df);
for i = 1:length(f_edge)-1
	sel		= (f >= f_edge(i)) & (f < f_edge(i+1));
	if sum(sel)
		avpwr_ft(i,:) = mean(pwr_ft(sel,:),1);
	end
end

% interpolate empty spaces in an intelligent manner (including [0,0])
empty = sum(avpwr_ft,2)==0;
if sum(empty)
	avpwr_ft(empty,:) = interp1([0; f_cntr(~empty)'], ...
		[zeros(1,size(avpwr_ft,2)); ...
		avpwr_ft(~empty,:)], f_cntr(empty)');
end
t = 1000*(t + t(2)/2);
f = f_cntr;

if nargout==0
	imagesc(t, f, 10*log10(avpwr_ft))
	axis xy
	colormap(jet)
	xlabel('Time (ms)')
	ylabel('Freq. (Hz)');
	fticks = linspace(f(1), f(end), length(f));
	n = ceil(length(f)/14);
	set(gca, 'TickDir','out', 'YTick', fticks(1:n:end), ...
		'YTickLabel', round(f(1:n:end)))
else
	logspec = avpwr_ft;
end

function [yo,fo,to] = specgram(varargin)
%SPECGRAM Spectrogram using a Short-Time Fourier Transform (STFT).
%   SPECGRAM has been replaced by SPECTROGRAM.  SPECGRAM still works but
%   may be removed in the future. Use SPECTROGRAM instead. Type help
%   SPECTROGRAM for details.
%
%   See also PERIODOGRAM, SPECTRUM/PERIODOGRAM, PWELCH, SPECTRUM/WELCH, GOERTZEL.

%   Author(s): L. Shure, 1-1-91
%              T. Krauss, 4-2-93, updated
%   Copyright 1988-2010 The MathWorks, Inc.
%   $Revision: 1.8.4.8 $  $Date: 2011/05/13 18:08:39 $

error(nargchk(1,5,nargin,'struct'));
[x,nfft,Fs,window,noverlap]=specgramchk(varargin);

nx = length(x);
nwind = length(window);
if nx < nwind    % zero-pad x if it has length less than the window length
	x(nwind)=0;  nx=nwind;
end
x = x(:); % make a column vector for ease later
window = window(:); % be consistent with data set

ncol = fix((nx-noverlap)/(nwind-noverlap));
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);
rowindex = (1:nwind)';
if length(x)<(nwind+colindex(ncol)-1)
	x(nwind+colindex(ncol)-1) = 0;   % zero-pad x
end

if length(nfft)>1
	df = diff(nfft);
	evenly_spaced = all(abs(df-df(1))/Fs<1e-12);  % evenly spaced flag (boolean)
	use_chirp = evenly_spaced & (length(nfft)>20);
else
	use_chirp = 0;
end

if (length(nfft)==1) || use_chirp
	y = zeros(nwind,ncol);
	
	% put x into columns of y with the proper offset
	% should be able to do this with fancy indexing!
	y(:) = x(rowindex(:,ones(1,ncol))+colindex(ones(nwind,1),:)-1);
	
	% Apply the window to the array of offset signal segments.
	y = window(:,ones(1,ncol)).*y;
	
	if ~use_chirp     % USE FFT
		% now fft y which does the columns
		y = fft(y,nfft);
		if ~any(any(imag(x)))    % x purely real
			if rem(nfft,2),    % nfft odd
				select = 1:(nfft+1)/2;
			else
				select = 1:nfft/2+1;
			end
			y = y(select,:);
		else
			select = 1:nfft;
		end
		f = (select - 1)'*Fs/nfft;
	else % USE CHIRP Z TRANSFORM
		f = nfft(:);
		f1 = f(1);
		f2 = f(end);
		m = length(f);
		w = exp(-1i*2*pi*(f2-f1)/(m*Fs));
		a = exp(1i*2*pi*f1/Fs);
		y = czt(y,m,w,a);
	end
else  % evaluate DFT on given set of frequencies
	f = nfft(:);
	q = nwind - noverlap;
	extras = floor(nwind/q);
	x = [zeros(q-rem(nwind,q)+1,1); x];
	% create windowed DTFT matrix (filter bank)
	D = window(:,ones(1,length(f))).*exp((-1i*2*pi/Fs*((nwind-1):-1:0)).'*f');
	y = upfirdn(x,D,1,q).';
	y(:,[1:extras+1 end-extras+1:end]) = [];
end

t = (colindex-1)'/Fs;

% take abs, and use image to display results
if nargout == 0
	newplot;
	if length(t)==1
		imagesc([0 1/f(2)],f,20*log10(abs(y)+eps));axis xy; colormap(jet)
	else
		% Shift time vector by half window length; the overlap factor has
		% already been accounted for in the colindex variable.
		t = ((colindex-1)+((nwind)/2)')/Fs;
		imagesc(t,f,20*log10(abs(y)+eps));axis xy; colormap(jet)
	end
	xlabel('Time')
	ylabel('Frequency')
elseif nargout == 1,
	yo = y;
elseif nargout == 2,
	yo = y;
	fo = f;
elseif nargout == 3,
	yo = y;
	fo = f;
	to = t;
end

function [x,nfft,Fs,window,noverlap] = specgramchk(P)
%SPECGRAMCHK Helper function for SPECGRAM.
%   SPECGRAMCHK(P) takes the cell array P and uses each cell as
%   an input argument.  Assumes P has between 1 and 5 elements.

x = P{1};
if (length(P) > 1) && ~isempty(P{2})
	nfft = P{2};
else
	nfft = min(length(x),256);
end
if (length(P) > 2) && ~isempty(P{3})
	Fs = P{3};
else
	Fs = 2;
end
if length(P) > 3 && ~isempty(P{4})
	window = P{4};
else
	if length(nfft) == 1
		window = hanning(nfft);
	else
		error(message('signal:specgram:NeedWindow'));
	end
end
if length(window) == 1, window = hanning(window); end
if (length(P) > 4) && ~isempty(P{5})
	noverlap = P{5};
else
	noverlap = ceil(length(window)/2);
end

% NOW do error checking
if (length(nfft)==1) && (nfft<length(window)),
	error(message('signal:specgram:WindowTooBig'));
end
if (noverlap >= length(window)),
	error(message('signal:specgram:OverlapTooBig'));
end
if (length(nfft)==1) && (nfft ~= abs(round(nfft)))
	error(message('signal:specgram:NFFTMustBePositive'));
end
if (noverlap ~= abs(round(noverlap))),
	error(message('signal:specgram:OverlapMustBePositive'));
end
if min(size(x))~=1,
	error(message('signal:specgram:MustBeVector'));
end

