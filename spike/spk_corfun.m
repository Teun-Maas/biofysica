function [XC,lags,ISI] = pa_spk_corfun(Spike,varargin)
% [XC,LAG,ISI] = PA_SPK_CORFUN(SPIKE)
%
% Correlation function.
%
% Obtain the autocorrelation function XC and the inter-spike-intervals ISI
% of structure SPIKE as obtained through PA_SPK_READSRC.
%
% PA_SPK_CORFUN(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'Fs'		- sample frequency (Hz). This is used to determine firing rate.
%	'display' - display the results
%
% See also PA_SPK_READSRC
%
% More information: "Spikes - exploring the neural code", Rieke et al. 1999, figure 2.5

% (c) 2012 Marc van Wanrooij
% E-mail: m.vanwanrooij@gmail.com

%% Variables
fac		= 10;
Ntrials = size(Spike,2);
Fs         = pa_keyval('Fs',varargin);
if isempty(Fs)
	Fs			= 1000; % (Hz)
end
%% Initialization
if isstruct(Spike) % if we have a structure containing Spike timings
	spikeon = pa_spk_timing_struct2mat(Spike,fac);
end


%% Plot variables
col         = pa_keyval('color',varargin);
if isempty(col)
	col			= 'k';
end
dspflag         = pa_keyval('display',varargin);
if isempty(dspflag)
	dspflag			= 0;
end


%% Do stuff
XC = 0;
for ii = 1:Ntrials
	%% Autocorrelation function
	t			= spikeon(ii,:);
	[C,lags]	= xcorr(t,t,50*fac,'coef');
	XC			= XC+C;
	
	%% Inter-spike interval
	t			= Spike(ii).spiketime;
	
% 	%% BrainWare check
% 	sw				= Spike(ii).spikewave;
% 	[~,mxindx]		= max(sw);
% 	[~,mnindx]		= min(sw);
% 	sel				= mxindx<mnindx;
% 	t(sel)			= t(sel)+0;
	
	
	t			= t/Fs*1000; % Spike times in (ms)
	t			= sort(t);
	isitmp		= diff(t);
	ISI(ii,1)	= {isitmp}; %#ok<AGROW>
end
lags	= lags/fac; % set lags to ms
XC		= XC/Ntrials*1000*fac; % normalize to rate


ISI		= horzcat(ISI{:});

if dspflag
	subplot(211)
	sel = XC==max(XC);
	mx = max(XC(~sel));
	plot(lags,smooth(XC),col);
	xlabel('\Delta t (ms)');
	ylabel('Conditional rate (spikes/s)');
	ylim([0 1.1*mx]);
	title('Autocorrelation');
	
	d	= 0.05;
	x	= -10:d:200;
	N	= hist(ISI,x);
	N	= N*1000/d/length(ISI);
	subplot(224)
	bar(x,N);
	xlim([0 50]);
	pa_verline;
	xlabel('\Delta t (ms)');
	ylabel('Conditional probability density (spikes/s)');
	title('Interval histogram');
end


