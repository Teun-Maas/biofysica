function Corr = pa_predictvocal(strfFname,vocFname,dname,sndFiles)

% PA_PREDICTVOCAL
%
% Assumptions: 0.25 octave, 2.5 octave bandwidth

%% Initialization
close all
if nargin<2
	vocFname = 'joe6711c01b00.mat';
	% 	vocFname = 'nep';
end
if nargin<1
	strfFname = 'joe6715c01b00.mat';
	% 	strfFname = 'nep';
end
if nargin<3
	dname = 'E:\DATA\Cortex\Test';
	dname2 = 'E:\DATA\Cortex\Test\Both';
end
if nargin<4
	sndFiles = 'E:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
	if ~exist(sndFiles,'file')
		sndFiles = 'C:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
	end
	if ~exist(sndFiles,'file')
		sndFiles = 'D:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
	end
	
end
msgid	= 'stats:nlinfit:IllConditionedJacobian';
warning('off',msgid);
dt		= 12.5; % (ms) STRF time bin = 1 / (velcoity step)
di		= 1;

% dt = 5;
% dspFlag = 1;

%% Load data
cd(dname);
vocFname = pa_fcheckexist(vocFname);
load(vocFname);

%% Spike density function ALL TRIALS
% Where does BF come from???

n = numel(spikeJ);
SND = NaN(n,1);
for ii = 1:n
	SND(ii) = spikeJ(ii).stimvalues(4);
end
uSND = unique(SND);
nSND = numel(uSND);
FR = struct([]);
FR2 = FR;
for ii = 1:nSND   % = 1:6, for all 6 sounds
	sel = SND == uSND(ii);
	s   = spikeJ(sel);
	nrep = sum(sel);
	MSDF	= pa_spk_sdf(s,'sigma',3,'Fs',1000,'kernel','gaussian');
	X		= 1:length(MSDF);
	% 	figure(666)
	% 	subplot(2,3,ii)
	% 	plot(X,MSDF,'k-')
	% 	hold on
	% 	xlabel('Time (ms)');
	% 	ylabel('Spike density / Firing rate (Hz)');
	% 	axis square;
	% 	box off
	% 	xlim([1 2500]);
	% 	ylim([0 400]);
	% 	pa_verline(300,'r-');
	% 	basemu = mean(MSDF(1:300));
	% 	pa_horline(basemu,'r-');
	
	MSDF		= MSDF(300:end);
	FR2(ii).sdf = MSDF;
	
	tim		= [s.spiketime];
	x		= 0:(dt/di):2500;
	fr		= hist(tim,x)/(dt/di)*1000/nrep;
	figure(666)
	subplot(2,3,ii)
	plot(x,fr,'k-');
	hold on
	xlabel('Time (ms)');
	ylabel('Spike density / Firing rate (Hz)');
	axis square;
	box off
	xlim([1 2500]);
	pa_verline(300,'r-');
	sel = x<300;
	basemu = mean(fr(sel));
	pa_horline(basemu,'r-');
	
	sel		= x>300;
	fr		= fr(sel);
	FR(ii).sdf = fr;
	
	
end

%% Reload data
cd(dname2);
strfFname = pa_fcheckexist(strfFname);
load(strfFname);

%% STRF
Pred = struct([]);
for ii = 1:6
	snd = ii;
	h = load(sndFiles);
	Y		= h.h{snd};
	dsp = 1;
	if dsp
		figure(668)
		subplot(3,2,ii)
		Fs = 50000;
		nsamples	= length(Y);
		wavplay(Y,Fs);
		t			= nsamples/Fs*1000;
		nseg		= t*di/dt;
		segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
		noverlap	= round(0.6*segsamples); % 1/3 overlap
		window		= segsamples+noverlap; % window size
		nfft		= 1000;
		xoct		= linspace(0.25,2.5,10);
		[S,F,T]		= spectrogram(Y,window,noverlap,nfft, Fs);
		P			= S.*conj(S);
		P			= log10(P);
		mxscal		= max((P(isfinite(P))));
		imagesc(T*1000,xoct,P);
		ylim([0.25 1.6]);
		set(gca,'YDir','normal');
		ylabel('Frequency (kHz)');
		caxis([-mxscal mxscal]);
		colorbar;
	end
	
	strf		= pa_spk_ripple2strf(spikeP,'comp',1);
	strf		= strf.strf;
	bf			= 1;
	if BF>250*(2^2.5)
		bf		= 2;
	end
	if BF>250*(2^5)
		bf		= 3;
	end
	bf			= 1;
	S			= pa_spk_predict(strf,Y,dt,'bf',bf,'display',ii,'interp',di,'abs',1);
	Pred(ii).S	= S;
end

%% Non-linearity
Pr = [];
R = [];
for ii = 1:6
	snd		= ii;
	
	a		= Pred(snd).S;
	n		= length(a);
	indx	= 1:n;
	b		= FR(snd).sdf;
	b		= b(indx);
	
	Pr		= [Pr a]; %#ok<AGROW>
	R		= [R b]; %#ok<AGROW>
end

figure(667)
plot(Pr,R,'k.')
hold on

beta0 = [500 -0.01 -150];
X = -600:600;
Y = pa_sigmoidfun(beta0,X);
plot(X,Y,'b-');

beta = nlinfit(Pr,R,@pa_sigmoidfun,beta0);

X = -600:600;
Y = pa_sigmoidfun(beta,X);
plot(X,Y,'r-');
axis square;

Corr = NaN(6,1);
for ii = 1:6
	snd = ii;
	S = Pred(ii).S;
	%% Equal lengths
	n = length(S);
	MSDF	= FR(snd).sdf;
	MSDF = MSDF(1:n);
	
	% 		S = pa_sigmoidfun(beta,S);
	figure(snd)
	subplot(326)
	cla;
	t = 1:length(S);
	t = t*12.5;
	plot(t,S,'k-','LineWidth',2);
	pa_horline;
	xlim([min(t) max(t)]);
	
	figure(snd)
	subplot(326)
	hold on
	X		= 1:length(MSDF);
	X		= X*12.5;
	plot(X,MSDF,'r-','LineWidth',2)
	
	
	%% Correlation
	[r,p] = corrcoef(S,MSDF(1:n));
	str = ['r = ' num2str(r(2),2) ', p =' num2str(p(2),2)];
	title(str);
	Corr(ii) = r(2);
	
	% 	subplot(326)
	% 	hold on
	% 	MSDF	= FR2(snd).sdf;
	% 	X		= 1:length(MSDF);
	% 	plot(X,MSDF,'g-','LineWidth',2)
	
	colorbar;
end

%%
warning('on',msgid);
