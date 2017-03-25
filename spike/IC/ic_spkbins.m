% Nspk=spkbins(S1,S2)
%   convert staircase spike-signal to spikes per sample
%

function Nspk = spkbins(S1,S2);

n = size(S1,1);
m = size(S1,2);

% combine two spike channels
x = zeros(2*m,n);
x(1:2:2*m-1,:) = S1';
x(2:2:2*m  ,:) = S2';
clear S1 S2

% calculate signal steps
Steps = [diff(x)' zeros(n,1)];
clear x

% calculate number of spikes in each bin
step = 65536/18;
Nspk = round(Steps/step);

% correct for reset jumps (16 levels)
s = Nspk<0;
Nspk(s) = Nspk(s)+17;


