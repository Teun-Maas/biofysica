function [HS,freq,theta] = modelHeadshadow(theta,freq,varargin)
% HS = MODELHEADSHADOW(THETA,FREQ)
%
% Approximation of the acoustic head shadow at horizontal angle THETA (deg)
% and frequency FREQ (Hz).
%
% By default this approximation is a single-pole, single zero filter, from:
% Brown CP, Duda RO. A structural model for binaural sound synthesis.
% IEEE Trans Speech Audio Process. 1998;6: 476–488. doi:10.1109/89.709673
%
% which works well for a spherical head.
%
% See also: modelILD

%% Initialization
% clearvars;
if nargin<2
	theta		= 0:5:180; % horizontal angle (deg)
% 	theta = [0 180];
	freq		= linspace(0,12000,1000); % frequency (Hz)
	[theta,freq] = meshgrid(theta,freq);
end

c			= keyval('c',varargin,343); % speed of sound (m/s)
a			= keyval('a',varargin,0.0875); % radius of head (m)
a			= keyval('a',varargin,0.57/(2*pi*2)); % radius of head (m)

minAlpha	= keyval('minAlpha',varargin,0.1); % minimum
minTheta	= keyval('minTheta',varargin,150); % minimum head shadow at minTheta, typically not at 180 deg due to summation of sounds coming from both sides
method		= keyval('method',varargin,'brown'); % use default Brown and Duda's filter approximation
method		= keyval('method',varargin,'rayleigh'); % use Rayleigh's head shadow physical model
method		= keyval('method',varargin,'duda'); % use Duda's description of Rayleigh's physical model
% method		= keyval('method',varargin,'cipic'); % use Duda's description of Rayleigh's physical model

threshold	= keyval('threshold',varargin,0.01); % threshold
r			= keyval('r',varargin,1); % distance (m)

graphFlag	= keyval('showPlot',varargin,false); % use default Brown and Duda's filter approximation
if nargin<1
graphFlag = true;
end

if exist('cbrewer','file')
cmap = cbrewer('div','RdYlBu',64,'pchip');
cmap = flipud(cmap);
else
	cmap = parula;
end
%% Reparametrization
o			= 2*pi*freq; % angular frequency (Hz)
o0			= c/a;
f0			= o0/(2*pi);

%% Head Shadow
switch method
	case 'brown'
		%% head shadow approximation with single-pole, single zero filter for spherical head
		% From: Brown CP, Duda RO. A structural model for binaural sound synthesis.
		% IEEE Trans Speech Audio Process. 1998;6: 476–488. doi:10.1109/89.709673
		mu			= 2*pi*a*freq/c;
		alpha		= (1+minAlpha/2)+(1-minAlpha/2).*cosd(theta./minTheta.*180);
		HS			= (1+1j.*alpha.*o./(2*o0))./(1+1j*o./(2*o0));
		HS			= 20*log10(abs(HS)); % head shadow in deciBe
	case 'opstal'
		%% fit sine through ILD data with square-root frequency dependency
		% from: van Opstal A.J. 2016, The Auditory System and Human
		% Sound-Localization Behavior, Elsevier, Academic Press 
		%
		% Note that azimuth/theta is defined differently than for Brown. 0
		% = straight ahead, +=right ear, -=left ear.
		HS				= -0.09*sqrt(freq).*sind(theta-90); 
	case 'cipic'
		[~,~,~,~,~,GI]	= getcipicavg;
		HS				= GI(freq,theta);
	case 'rayleigh'
		HS = getrayleigh;
		return
	case 'duda'
		 HS = sphere(theta, freq,'c',c,'a',a,'r',r);
end

%% Graphics
if graphFlag
close all
colormap(cmap);
	% figure(1)
	% clf
	% subplot(221)
	% plot(theta',alpha','ko-',...
	% 	'MarkerFaceColor','w');
	% nicegraph;
	
	% % xi = oct2bw(1000,-1:4);
	% xi = [0.1 1 10 100];
	% subplot(122)
	% semilogx(mu,HS,'k-',...
	% 	'MarkerFaceColor','w');
	% set(gca,'XTick',xi,'XTickLabel',xi)
	% nicegraph;
	% ylim([-20 10]);
	% xlim([0.1 100]);
	subplot(131)
	semilogx(freq,HS,'-',...
		'MarkerFaceColor','w','LineWidth',2);
	xi	= oct2bw(1000,-2:4);
	set(gca,'XTick',xi,'XTickLabel',xi/1000)
	nicegraph;
	ylim([-20 10]);
	xlim([100 12000]);
	verline(f0);
	ylabel('response (dB)');
	xlabel('frequency (kHz)');
	
	subplot(132)
	plot(theta',HS','-',...
		'MarkerFaceColor','w','LineWidth',2);
	nicegraph;
	ylim([-20 10]);
	xlim([0 180]);
	ylabel('response (dB)');
	xlabel('azimuth (deg)');
	xi	= 0:30:180;
	set(gca,'XTick',xi,'XTickLabel',xi)
	
	
	x = freq(:,1)';
	y = theta(1,:)'-90;
	z = HS';
		subplot(133)
		imagesc(x,y,z);
		hold on
		[C,h] = contour(x,y,z,20,'w');

		nicegraph;
% 	plot(theta',HS','-',...
% 		'MarkerFaceColor','w','LineWidth',2);
% 	nicegraph;
% 	ylim([-20 10]);
% 	xlim([0 180]);
% 	ylabel('response (dB)');
% 	xlabel('azimuth (deg)');
% 	xi	= 0:30:180;
% 	set(gca,'XTick',xi,'XTickLabel',xi)
	% subplot(224)
	% contourf(freq(:,1),theta(1,:),HS',10);
	% nicegraph;
end

if nargout<1
	clearvars
end


function [HS_R,HS_L,ILD,freq,theta,F] = getcipicavg

if ~exist('cipicavg.mat','file')
	load('cipicavg.mat','HS_L','HS_R','ILD','freq','theta','F');
else
	
	fname = which('anthro.mat');
	load(fname,'id');
	
	M			= zeros(5,24);
	HS			= M;
	nid			= length(id);
	MILD		= 0;
	MAR			= 0;
	MAL			= 0;
	start_flag	= 0;             % Prevents any action until data is available
	log_scale	= 0;             % Default to use a logarithmic frequency scale
	do_smooth	= 1;             % Default to use constant-Q smoothing
	Q			= 8;             % For constant-Q smoothing, of course
	fs			= 44100;         % Sampling frequency in Hz
	fNy			= fs/2000;       % Nyquist frequency in kHz
	Ia			= 13;            % Initial azimuth index   (0 deg)
	Ie			= 9;             % Initial elevation index (0 deg)
	dBmin		= -30;                % Lower limit for spectral plots in dB
	dBmax		= +20;                % Upper limit for spectral plots in dB
	Fmin = 0.1; Fmax = 20;      % Min and max frequency in kHz
	Tmin = 0.3; Tmax = 2.5;     % Min and max time in ms
	LT			= 200;                                            % Number of time points
	T			= (0:(1/fs):((LT-1)/fs))*1000;                        % Times in ms
	miI			=  0.0; mxI = 0.75;     % Min and max ITD in ms
	miT			= find(T<=Tmin,1,'last');
	mxT			= find(T>=Tmax,1,'first');
	NFFT		= LT;                  % More freq resolution if NFFT > LT, but slows computing
	
	A			= [-80 -65 -55 -45:5:45 55 65 80];   LA = length(A);  % Azimuths
	E			= -45:(360/64):235;                  LE = length(E);  % Elevation
	for ii = 1:nid
		sub			= id(ii);
		sel			= id==sub;
		subject_numbers = sub_nums; % Load cell array of subject numbers from file
		
		
		%% Data
		c = ['CIPIC/standard_hrir_database/subject_' num2str(sub,'%03d')];
		d = what(c);
		c = d.path;
		cd(c);
		load('hrir_final','hrir_l','hrir_r');
		hl		= squeeze(hrir_l(:,Ie,:))';               % Normalize data [-1 1]
		hr		= squeeze(hrir_r(:,Ie,:))';
		[AL]	= freq_resp(hl,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
		[AR,F] = freq_resp(hr,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
		ILD		= AR'-AL';
		MILD	= MILD+ILD;
		MAL		= MAL+AL;
		MAR		= MAR+AR;
	end
	
	ILD		= MILD/nid;
	HS_R	= MAR/nid;
	HS_R	= HS_R-mean(HS_R,2); % directional component only
	HS_L	= MAL/nid;
	HS_L	= HS_L-mean(HS_L,2); % directional component only
	
	freq	= F;
	theta	= -(A+90)+180;
	theta = fliplr(theta);
	HS_R = fliplr(HS_R);
	
	[freq,theta] = ndgrid(freq,theta);
	
	F = griddedInterpolant(freq,theta,HS_R,'linear');
	
	c = 'biofysica/audition';
	d = what(c);
	c = d.path;
	cd(c);
	save('cipicavg','HS_R','HS_L','ILD','freq','theta','F');
end


function L = getrayleigh
%% Work In Progress

%% General
% lamba = wavelength
% c		= radius sphere
% k		= 2pi/lambda
% mu	= cos(theta)
% I = intensity = F^2+G^2

kc = [0.5 1 2];
mu = [1 -1 0];
[kc,mu] = meshgrid(kc,mu);

I = [0.294291 0.259729 0.231999;
	0.502961 0.285220 0.236828;
	0.6898 0.3182 0.3562]';

figure(1)
clf
plot(kc',I','o-',...
'MarkerFaceColor','w','LineWidth',2,'MarkerSize',15)
nicegraph;
xlabel('kc');
ylabel('intensity');
xlim([0 2.5]);
legend(num2str(mu(:,1)));

%% kc = 2, more azimuth angles
% Note that the sound source was not softest at the anti-pole
% (mu=-1,theta=180 deg), but some 40 deg off from this angle.
% 
% Determined with Legendre series, up to P6
theta = 0:15:180;
FG = [
	0.7968+0.2342i
	0.8021+0.1775i
	0.7922+0.0147i
	0.7139-0.2287i
	0.5114-0.4793i
	0.1898-0.6247i
	-0.1538-0.5766i
	-0.3790-0.3413i
	-0.3992-0.0243i
	-0.2401+0.2489i
	-0.0088+0.4157i
	0.1781+0.4883i
	0.2495+0.5059i];

M = abs(FG);
I = M.^2;
I = 4*I; % for large spheres, the sphere acts as a plane, so intensity is quadrupled
L = 20*log10(I);

figure(2)
clf
plot(theta,L,'ko-',...
'MarkerFaceColor','w','LineWidth',2,'MarkerSize',15)
nicegraph;
xlabel('azimuth \theta (deg)');
ylabel('level (dB)');
set(gca,'XTick',0:45:180);
xlim([-30 210]);



function HS = sphere(Theta, Freq,varargin)
% HS = SPHERE(A,R,THETA,F,C,THSRESHSOLD)
%
% Head shadow for sphere with radius A.
% c		= ambient speed of sound
% freq	= frequency (Hz)
% theta = angle of incidence (deg)
% r		= distance from the center of the sphere to the source (m)
%
% Iteration stops when the fractional change falls below the THRESHOLD for
% two successive terms.  
%
% From: Duda RO, Martens WL. Range dependence of the response of a spherical
% head model. J Acoust Soc Am. Acoustical Society of America; 1998;104:
% 3048–3058. doi:10.1121/1.423886  


%% Initialization
if nargin<2
	theta		= 0:15:180; % horizontal angle (deg)
	freq		= linspace(0,12000,1000); % frequency (Hz)
	[Theta,Freq] = meshgrid(theta,freq);
end
c			= keyval('c',varargin,343); % speed of sound (m/s)
a			= keyval('a',varargin,0.0875); % radius of head (m)
threshold	= keyval('threshold',varargin,0.01); % threshold
r			= keyval('r',varargin,10); % distance (m)

[m,n]		= size(Theta);


%%
HS			= NaN(m,n);
for ff = 1:m
	for tt = 1:n
		theta		= Theta(ff,tt);
		freq		= Freq(ff,tt);
		x			= cosd(theta);
		mu			= (2.*pi.*freq.*a)./c;
		rho			= r./a;
		zr			= 1./(1i.*mu.*rho);
		za			= 1./(1i.*mu);
		Qr2			= zr;
		Qr1			= zr.*(1-zr);
		Qa2			= za;
		Qa1			= za.*(1-za);
		P2			= 1;
		P1			= x;
		s			= 0;
		term		= zr./(za.*(za - 1));
		s			= s+term;
		term		= (3.*x.*zr.*(zr-1))./(za.*(2.*za.^2-2.*za+1));
		s			= s+term;
		oldratio	= 1;
		newratio	= abs(term)./abs(s);
		m			= 2;
		while (oldratio > threshold) || (newratio > threshold)
			Qr		= -(2.*m-1).*zr.*Qr1+Qr2;
			Qa		= -(2.*m-1).*za.* Qa1+Qa2;
			P		= ((2.*m-1).*x.*P1-(m-1).*P2)./m;
			term	= ((2.* m+1).*P.*Qr)./((m+1).*za.*Qa-Qa1);
			s		= s+term;
			m		= m+1;
			Qr2		= Qr1;
			Qr1		= Qr;
			Qa2		= Qa1;
			Qa1		= Qa;
			P2		= P1;
			P1		= P;
			oldratio = newratio;
			newratio = abs(term)./abs(s);
		end
		HS(ff,tt)			= (rho.*exp(-1i*mu).*s)./(1i.*mu);
	end
end

%%
HS = 20*log10(abs(HS));

