function pa_spk_plotstrf(data)
% PA_SPK_PLOTSTRF(STRF)
%
% Plot Spectro-Temporal Receptive Field in structure STRF, ontained by
% PA_SPK_RIPPLESTRF. This structure contains the following fields:
% - magnitude
% - phase
% - density
% - velocity
% - norm
% - strf
% Additional fields (depending on experiment) are:
% - magextra
% - phasextra
%
%
% See also PA_SPK_RIPPLE2STRF

% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Initialization
% Data-structure contains N elements, one for each Modulation Depth
% presented.
n = length(data);

%% Graphics
% Plot magnitude and phase plots
for ii = 1:n
	figure(ii)
	m		= data(ii).magnitude;
	p		= data(ii).phase;
	nrm		= data(ii).norm;
	uD		= data(ii).density;
	uV		= data(ii).velocity;
	plotMP(m,p,nrm,uD,uV);
end

% Sometimes only a subsection of ripples were presented, with a fixed
% density and velocity, instead of an entire quadrant.
if isfield(data,'magextra')
	for ii = 1:n
		figure(100+ii)
		m		= data(ii).magextra;
		p		= data(ii).phasextra;
		nrm		= data(ii).norm;
		uD		= data(ii).density;
		uV		= data(ii).velocity;
		plotMP(m,p,nrm,uD,uV);
	end
end



%% STRF
	nsb = ceil(sqrt(n)); %number of subplots
for ii = 1:n
	figure(200)
	subplot(nsb,nsb,ii)
	x		= data(ii).time;
	dx = mean(diff(x));
	y		= data(ii).frequency;
	z		= data(ii).strf;
	[~,m] = max(max(z,[],1));
	
	mxscal = max(abs(z(:)));
	XI = linspace(min(x),max(x),500);
	YI = linspace(min(y),max(y),500);
	ZI = interp2(x,y',z,XI,YI','*cubic');
	
	subplot(121)
	imagesc(x,y,z);
	colorbar;
	axis square
	caxis([-mxscal,mxscal]);
	ylabel('Frequency (octaves)');
	xlabel('Time (ms)');
	title('STRF');
	set(gca,'YDir','normal','TickDir','out');
	pa_verline(x(m));
	xlim([0 100])
	
	subplot(122)
	imagesc(XI,YI,ZI);
	shading flat;
	colorbar;
	axis square
	caxis([-mxscal,mxscal]);
	ylabel('Frequency (octaves)');
	xlabel('Time (ms)');
	title('STRF');
	set(gca,'YDir','normal','TickDir','out');
	pa_verline(x(m));
	xlim([0 100])

% 	col = pa_statcolor(64,[],[],[],'def',8,'disp',false);
% 	colormap(col);
end

function plotMP(M,P,Nrm,x,y)
% Magnitude plot
subplot(131)
z		= M;
imagesc(x,y,z);
colorbar;
ylabel('ripple velocity (Hz)');
axis square
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',y,'XTick',x);
title('Magnitude (Hz)');


% Phase plot
subplot(132)
z		= P;
imagesc(x,y,z);
caxis([-pi pi])
colorbar;
xlabel('ripple density (cycles/s)');
axis square
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',[],'XTick',x);
title('Phase');

% Norm plot
subplot(133)
z		= Nrm;
imagesc(x,y,z);
caxis([0 1])
colorbar;
axis square
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',[],'XTick',x);
title('Norm factor');
% try
% 	col = pa_statcolor(64,[],[],[],'def',7,'disp',false);
% 	colormap(col);
% 	end