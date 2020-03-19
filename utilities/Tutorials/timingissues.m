
% see also Matlab for neuroscientists Chapter 6
close all
clearvars;
format long; % better time resolution
% format short; % default
n = 100000;
%% no preallocation


ii = 0;
t = [];
while ii<n+1
	tic;
	ii = ii+1;
% 	dur = 0.4+0.05*randn(1);
	dur = 0;
% 	pause(dur);
	t(ii,1) = toc-dur;
end
figure(1)
clf
ax(1) = subplot(211);
plot(t(50:end)*1000)
title(mean(t(50:end)*1000));

%% with preallocation
ii = 0;
t = NaN(n,1);

while ii<n+1
	tic;
	ii = ii+1;
% 	dur = 0.4+0.05*randn(1);
	dur = 0;
% 	pause(dur);
	t(ii,1) = toc-dur;
end
figure(1)
ax(2) = subplot(212);
plot(t(50:end)*1000)
title(mean(t(50:end)*1000));

linkaxes(ax);
xlabel('Trial #');
ylabel('Time (ms)');
max(t*1000)