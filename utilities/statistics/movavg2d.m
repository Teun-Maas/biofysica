function [avg] = movavg2d(data,windowSize)
% [AVG, XAVG] = MOVAVG(DATA,WINDOWSIZE);
%
% Determines a moving average, AVG, of continuous DATA
% with a window of WINDOWSIZE samples (default = 20).
%
% See also FILTER



%% Initialization

if nargin<2
    windowSize = min(20,size(data,1));
end
if nargin<3
    x   = 1:length(data);
end

M = data;
for ii = 1:size(M,1)
	a = M(ii,:);
	b = movavg(a,windowSize);
	M(ii,:) = b;
end

for ii = 1:size(M,2)
	a = M(:,ii);
	b = movavg(a,windowSize);
	M(:,ii) = b;
end
avg = M;
