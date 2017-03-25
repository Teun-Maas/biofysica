function spikemat = pa_spk_timing_struct2mat(Spike,fac)
% SPIKEMAT = PA_SPK_TIMING_STRUCT2MAT(SPIKE)
%
% Convert a strcuture SPIKE containing spike timings in field 'spiketime'
% to a matrix SPIKEMAT of ones and zeros, with ones indicating a spike
% present at a time denoted by the matrix index.
%
% SPIKEMAT = PA_SPK_TIMING_STRUCT2MAT(SPIKE,FAC)
%
% Additonal input argument FAC 'increases' the sampling rate.
% By default the matrix's sampling rate is defined by rounding the spike
% timings to integers (for BrainWare this defaults to a ms-range, FAC=1).
%
% See also PA_SPK_SDF, PA_SPK_CORFUN

% 2012 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<2
	fac = 1;
end

% Spike timing is set in ms
% Increase sampling rate by a factor
Ntrials		= size(Spike,2);
spikemat		= zeros(Ntrials,ceil(max([Spike.spiketime]*fac)));
for ii = 1:Ntrials
	timings			= Spike(ii).spiketime;
	
% 	%% BrainWare check
% 	spikes			= Spike(ii).spikewave;
% 	[~,mxindx]		= max(spikes);
% 	[~,mnindx]		= min(spikes);
% 	sel				= mxindx<mnindx;
% 	timings(sel)	= timings(sel)+1;
	
	indx			= ceil(timings*fac);
	spikemat(ii,indx) = 1;
end
