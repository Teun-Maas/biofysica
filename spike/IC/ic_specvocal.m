function s = ic_specvocal(x,f1,f2,df,dt)
% S = IC_SPECVOCAL(x,<f1>,<f2>,<df>,<dt>);
%
%	s{1..6}: structure with 6 spectra
%
%	x: structure with timetraces of calls,
%		e.g. load timetraces
%	f1, f2 lower and upper boundary of spectral window
%	df, dt: spectral and temporal bin
%	defaults: f1=297 Hz, f2=1420 Hz, df=0.25 oct, dt=12.5 ms
%	alternative: shift of spectrum, uses specgram3
%

% Huib Versnel/John van Opstal
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<5
	dt = 12.5; % (ms) = 1/(velocity step)
end
if nargin<4
	df = 0.25; % (oct) = 1/(density step)
end
if nargin<3
	f2 = 250*2^2.5; % 1414 Hz
end
if nargin<2
	f1 = 250*2^0.25; % 297 Hz
end
Nf			= round(log2(f2/f1)/df)+1;
freq		= (0:Nf-1)*df;

s	= cell(6,1);
for ii = 1:6,
	t		= x{ii};
	f		= ic_specgram(t,[f1 f2],df,dt);
	P		= 10*log10(f);
	s{ii}	= f;
	Nt		= size(f,2);
	t		= (0:Nt-1)*dt;
	
	%% Graphics
	subplot(3,2,ii);
	imagesc(t,freq,P);
	set(gca,'YDir','normal','TickDir','out');
	colorbar;
% 	axis square
end;
