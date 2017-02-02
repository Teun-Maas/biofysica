function bw = pa_freq2bw(F1,F2)
% BW = PA_FREQ2BW(F1,F2)
%
% Determine bandwidth in octaves between low frequency F1 and high
% frequency F2


% (c) 2012 Marc van Wanrooij

bw = F2/F1; % bandwidth in Hz
bw = log2(bw); % bandwidth in octaves
