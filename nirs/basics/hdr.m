function Y = hdr(X,varargin)
% Y = HDR(BETA,X)
%
% Canonical hemodynamic response function (SPMs double gamma function)
%
% See also:
% - https://en.wikipedia.org/wiki/Haemodynamic_response
% - Lindquist, M. A., Meng Loh, J., Atlas, L. Y., & Wager, T. D. (2009). 
% Modeling the hemodynamic response function in fMRI: efficiency, bias and
% mis-modeling. NeuroImage, 45(1 Suppl), S187-98.
% http://doi.org/10.1016/j.neuroimage.2008.10.065  

%% Initialization
Fs = keyval('Fs',varargin,250); % sample frequency (Hz)/(samples/s)
if nargin<1
	X		= zeros(50*Fs,1);
	X(round(3*Fs):round(23*Fs)) = 1;
end
dur			= keyval('dur',varargin,50); % duration (s)
dursamples	= round(dur*Fs);
dispFlag    = keyval('disp',varargin,true); % duration (s)

%% Hemodynamic impulse response


N       = length(X); % samples
t		= 0:dursamples-1; % the single hemodynamic response time (samples)
t		= t/Fs; % s
a1		= 6; % peak time 4.5 s after stim, shape
a2		= 16; % peak time 4.5 s after stim, shape
b1		= 1;
b2		= 1;
c		= 1/6;
y1       = gampdf(t,a1,b1);
y2		= c*gampdf(t,a2,b2);
y		= y1-y2;
Y       = conv(X,y);
Y       = Y(1:N);

t = 0:N-1; % time (for entire trace)
t = t/Fs;
if dispFlag



end

if dispFlag
	graph(t,X,Y,N,Fs);
end

function graph(t,X,Y,N,Fs)

Y = Y/max(Y);
	figure(1)
	clf
	plot(t-3,X,'k-','LineWidth',2);
	
	hold on
	plot(t-3,Y,'LineWidth',2);
	plot(t(1)-3,Y(1),'ko','MarkerFaceColor','w');
	plot(t(end)-3,Y(end),'ko','MarkerFaceColor','w');
	
	axis square;
	ylim([-0.2 1.2]);
	xlabel('Time (s)');
	ylabel('Amplitude (au)');
	set(gca,'TickDir','out');
	box off
	verline([0 5]);
	verline([1 2 3 4 12],'k:');
	horline([0:0.1:1],'k:');
% 	xlim([0 5])
set(gca,'XTick',(0:1:N/Fs)-3);


	try
	nicegraph;
	savegraph(mfilename,'eps');
	catch
		disp('Why don''t you try out the Biofysica Toolbox @ Gitlab.ru.nl');
	end