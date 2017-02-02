function  dSPL = pa_addspl(varargin)
% DSPL = PA_ADDSPL
%
% Determine SPL difference when two sounds are added.

grph         = pa_keyval('display',varargin);
if isempty(grph)
	grph			= 1;
end

%% Tones
Beta1	= 0;
Beta2	= 0;
f1		= 100; % (Hz)
f2		= 200;
a1		= 1;
a2		= 1;
stm1		= pa_gentone(f1,100,'display',0,'phase',pa_deg2rad(Beta1),'nenvelope',0);
[stm2,Fs]	= pa_gentone(f2,100,'display',0,'phase',pa_deg2rad(Beta2),'nenvelope',0);
stm1		= a1*stm1;
stm2		= a2*stm2;

%% Derivatives
df		= abs(f2-f1); % Hz
N		= length(stm1); % samples
t		= (1:N)/Fs; % (s)

%% Parameters
pSnd1	= rms(stm1); % average pressure sound 1
pSnd2	= rms(stm2);
pSnd12	= rms(stm1)*rms(stm2);
pMax	= max([pSnd1 pSnd2]);

if df == 0
	%% Addition of coherent sound pressures (Bies and Hansen, Ch. 1, 1996)
	% if two sounds of the same frequency are added, the phase difference
	% between the two sounds has to be included
	% (CH Hansen, Fundamentals of Acoustics: http://www.who.int/occupational_health/publications/noise1.pdf)
	Dphase	= abs(pa_deg2rad(Beta2)-pa_deg2rad(Beta1));
	p2Total	=  pSnd1^2 + pSnd2^2 + 2*(pSnd12)*cos(Dphase);
	dSPL	= 10*log10(p2Total/pMax^2);
	
else
	%% Beating
	p2Total	= pSnd1^2 + pSnd2^2 + 2*(pSnd12)*cos(2*pi*df*t); % dynamically changing
	pbeat	= sqrt(2)*sqrt(p2Total);
	p2Total	= rms(p2Total);
	dSPL	= 10*log10(p2Total/pMax^2)
	
	%% Addition of incoherent sound pressures (logarithmic addition)
	p2Total	= pSnd1^2 + pSnd2^2; % limit
	dSPL	= 10*log10(p2Total/pMax^2)
end

if grph
	t = t*1000;
	close all
	subplot(311)
	plot(t,stm1)
	hold on
	h = pa_horline(sqrt(2)*pSnd1,'k-');
	set(h,'LineWidth',2);
	h = pa_horline(-sqrt(2)*pSnd1,'k-');
	set(h,'LineWidth',2);

	subplot(312)
	plot(t,stm2)
	hold on
	h = pa_horline(sqrt(2)*pSnd2,'k-');
	set(h,'LineWidth',2);
	h = pa_horline(-sqrt(2)*pSnd2,'k-');
	set(h,'LineWidth',2);
	subplot(313)
	plot(t,stm1+stm2)
	hold on
	if exist('pbeat','var')
		plot(t,pbeat,'k-','LineWidth',2);
		plot(t,-pbeat,'k-','LineWidth',2);
	else
		h = pa_horline(sqrt(2)*sqrt(p2Total),'k-');
		set(h,'LineWidth',2);
		h = pa_horline(-sqrt(2)*sqrt(p2Total),'k-');
		set(h,'LineWidth',2);
		
	end
	xlabel('Time (ms)');
end
