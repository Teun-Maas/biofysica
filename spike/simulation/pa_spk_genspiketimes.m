function T = pa_genspiketimes(N, duration)
% TIM = PA_GENSPIKETIMES(N,DUR)
%
% Generate random spike Tings for N spikes within a signal of duration
% DUR
%
% See also PA_SPIKEWAVE

% (c) 2011 Marc van Wanrooij

T		= floor(duration*rand(N,1))+1;
T		= sort(T);
