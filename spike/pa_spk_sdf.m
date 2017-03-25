function [MSDF,sdf] = pa_spk_sdf(spike,varargin)
% [MSDF,SDF] = PA_SPK_SDF(SPIKE)
%
% Obtain (mean) spike-density function (M)SDF by convolving trains of action
% potentials in SPIKE with a Gaussian, Possion of growth-decay function
% kernel.
%
% SPIKE should be a (sparse) MxN matrix of M trials and N samples,
% consisting of ones (spike) and zeros (no spikes). SPIKE can also be a
% structure containing spike timings in the field spiketime. MSDF will be 
% the average across trials, while SDF will also be a MxN matrix.
%
% PA_SDF(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'kernel'	- specify convolution kernel. Options are:
%			- 'gaussian' (default)
%			- 'poisson'
%			- 'EPSP' a combination of growth and decay funcions that
%			resemble a postsynaptic potential. Default growth and decay
%			time constant: 1 and 20 (chosen according to Thompson et al. 1996)
%
% And some kernel parameters might be requested:
%	'sigma'	- specify Gaussian kernel width. Default: 5 samples.
%	'Te'	and 'Td' - specify growth and decay time constants for EPSP
%						kernel. Default: 1 and 20.
%
% To obtain a firing rate instead of a density function, input a sampling
% frequency:
%	'Fs'	- sample frequency (Hz). This is used to determine firing rate.
%			By default Fs will be empty and only the spike density function
%			is determined.
%
% See also PA_SPK_RASTERPLOT
%
% More information: 
%		"Spikes - exploring the neural code", Rieke et al. 1999, figure	2.1
%		"Matlab for Neuroscientists", Wallisch et al. 2009,  section 15.5.2


% (c) 2011 Marc van Wanrooij
% E-mail: m.vanwanrooij@gmail.com

%% Initialization
if isstruct(spike) % if we have a structure containing spike timings
	spike = pa_spk_timing_struct2mat(spike);
end

if issparse(spike) % if we have a sparse matrix
	spike = full(spike); % make it full
end

%% Optional input
sigma         = pa_keyval('sigma',varargin);
if isempty(sigma)
	sigma			= 5;
end
Te         = pa_keyval('Te',varargin);
if isempty(Te)
	Te			= 1;
end
Td         = pa_keyval('Td',varargin);
if isempty(Td)
	Td			= 20;
end
kernel         = pa_keyval('kernel',varargin);
if isempty(kernel)
	kernel			= 'gaussian';
end
Fs         = pa_keyval('Fs',varargin);

[ntrials,nsamples]		= size(spike);
% if ntrials>nsamples
% 	disp('More trials than samples');
% 	disp('Is this correct?');
% end
switch kernel
	case 'gaussian'
		winsize		= sigma*5;
		t			= -winsize:winsize;
		window		= normpdf(t,0,sigma);
		winsize		= [winsize nsamples+winsize-1];
	case 'poisson'
		winsize		= [1 nsamples];
		t			= 1:150;
		window		= poisspdf(t,sigma);
	case 'EPSP'
		winsize		= [1 nsamples];
		t			= 0:Td*10;
		window		= (1-exp(-t/Te)).*(exp(-t/Td));
		window		= window./sum(window);
	case 'pulse'
		winsize		= [1 nsamples];
		window = pa_spk_spikewave(1);
	case 'boxcar'
		winsize		= sigma;
		t			= -winsize:winsize;
		window		= ones(size(t));
		winsize		= [winsize nsamples+winsize-1];		
	otherwise % also do 'gaussian'
		winsize		= sigma*5;
		window		= normpdf(-winsize:winsize,0,sigma);
		winsize = [winsize nsamples+winsize-1];
end


%% Convolute
sdf			= spike;
for ii = 1:ntrials
	convspike	= conv(spike(ii,:),window);
	convspike	= convspike(winsize(1):winsize(2));
	sdf(ii,:)	= convspike;
end

%% Obtain firing rate
if ~isempty(Fs) % if sampling frequency is given
	MSDF		= sum(sdf)/ntrials*Fs; % Firing Rate
	sdf			= sdf*Fs;
	if strcmp(kernel,'boxcar')
		sdf = sdf/(sigma*2);
	end
else % just average
	MSDF		= sum(sdf)/ntrials; % Probability
end