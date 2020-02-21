function ILD = modelILD(theta,freq,varargin)

% ILD = MODELILD(THETA,FREQ)
%
% Approximation of the ILD at horizontal angle THETA (deg)
% and frequency FREQ (Hz).
%
% By default this approximation is from a single-pole, single zero head
% shadow filter, from: 
% Brown CP, Duda RO. A structural model for binaural sound synthesis.
% IEEE Trans Speech Audio Process. 1998;6: 476–488. doi:10.1109/89.709673
% which works well for a spherical head.
%
% See also: MODELHEADSHADOW

%% Initialization
% clearvars;
if nargin<2
	theta		= 0:15:180; % horizontal angle (deg)
	freq		= linspace(0,16000,1000); % frequency (Hz)
	[theta,freq] = meshgrid(theta,freq);
end

c			= keyval('c',varargin,343); % speed of sound (m/s)
a			= keyval('a',varargin,0.0875); % radius of head (m)
minAlpha	= keyval('minAlpha',varargin,0.1); % minimum
minTheta	= keyval('minTheta',varargin,150); % minimum head shadow at minTheta, typically not at 180 deg due to summation of sounds coming from both sides
method		= keyval('method',varargin,'brown'); % use default Brown and Duda's filter approximation
% method		= keyval('method',varargin,'opstal'); % use default Brown and Duda's filter approximation
method		= keyval('method',varargin,'cipic'); % use default Brown and Duda's filter approximation
method		= keyval('method',varargin,'duda'); % use default Brown and Duda's filter approximation
graphFlag	= keyval('showPlot',varargin,false); % plot
if nargin<1
graphFlag = true;
end
%% Reparametrization
o			= 2*pi*freq; % angular frequency (Hz)
o0			= c/a;
f0			= o0/(2*pi);

%% Head Shadow
% switch method
% 	case 'brown'
		%% head shadow approximation with single-pole, single zero filter for spherical head
		% From: Brown CP, Duda RO. A structural model for binaural sound synthesis.
		% IEEE Trans Speech Audio Process. 1998;6: 476–488. doi:10.1109/89.709673
		HS_R = modelHeadshadow(theta,freq,...
			'c',c,'a',a,'minAlpha',minAlpha,'minTheta',minTheta,'method',method,'showPlot',false);

		HS_L = modelHeadshadow(180-theta,freq,...
			'c',c,'a',a,'minAlpha',minAlpha,'minTheta',minTheta,'method',method,'showPlot',false);
% end


%% ILD
% A positive ILD corresponds to right ear loudest
ILD = HS_R-HS_L;

%% Graphics
if graphFlag
	figure(1)
	clf
	subplot(121)
	semilogx(freq,ILD,'-',...
		'MarkerFaceColor','w','LineWidth',2);
	xi	= oct2bw(1000,-2:4);
	set(gca,'XTick',xi,'XTickLabel',xi/1000)
	nicegraph;
	ylim([-20 20]);
	xlim([100 12000]);
	verline(f0);
	ylabel('ILD (dB)');
	xlabel('frequency (Hz)');
	
	subplot(122)
	plot(theta(1:100:end,:)',ILD(1:100:end,:)','o-',...
		'MarkerFaceColor','w','LineWidth',2);
	nicegraph;
	ylim([-20 20]);
	xlim([0 180]);
	ylabel('ILD (dB)');
	xlabel('azimuth (deg)');	
	xi	= 0:30:180;
	set(gca,'XTick',xi,'XTickLabel',xi)
	
	figure(2)
	clf
	subplot(131)
	contourf(freq,theta,HS_R,20);
	nicegraph;
	xi	= oct2bw(1000,-2:4);
	set(gca,'XTick',xi,'XTickLabel',xi/1000,'XScale','log');
	caxis([-10 10]);

	subplot(132)
	contourf(freq,theta,HS_L,20);
	nicegraph;
	xi	= oct2bw(1000,-2:4);
	set(gca,'XTick',xi,'XTickLabel',xi/1000,'XScale','log');
	caxis([-10 10]);
	
	subplot(133)
	contourf(freq,theta,ILD,20);
	nicegraph;
	xi	= oct2bw(1000,-2:4);
	set(gca,'XTick',xi,'XTickLabel',xi/1000,'XScale','log');
	caxis([-10 10]);
end


if nargout<1
	clearvars
end