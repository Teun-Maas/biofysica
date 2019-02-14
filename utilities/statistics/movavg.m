function [avg,xavg] = movavg(data,windowSize,x)
% [AVG, XAVG] = MOVAVG(DATA,WINDOWSIZE);
%
% Determines a moving average, AVG, of continuous DATA
% with a window of WINDOWSIZE samples (default = 20).
%
% See also FILTER



%% Initialization
if nargin<1
	% just some test
    data    = normpdf(-100:100,30,30);
    noise   = randn(1,length(data))/1000;
    data    = data+noise;
	% try this function by using pa_runavg without in- or output
end
if nargin<2
    windowSize = 20;
end
if nargin<3
    x   = 1:length(data);
end

%% Find Running Average by Filter-method
% This function uses FILTER instead of a FOR-loop,
% and thus should be a lot faster.
avg     = filtfilt(ones(1,windowSize)/windowSize,1,data);

%% Correct the independent variable
xavg    = x-(windowSize./2).*mean(diff(x));

%% Some graphics if there is no other output
if nargout<1
    close all;
    h1 = plot(x,data,'ko'); set(h1,'MarkerFaceColor','w');
    hold on;
    h2 = plot(xavg,avg,'r-'); set(h2,'LineWidth',3);
    xlabel('Sample Number');
    ylabel('Amplitude');
end