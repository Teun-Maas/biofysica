function [M,S] = pa_geninfomask(Par)
% [MASKER,TARGET] = PA_GENINFOMASK
%
% Generate MASKER and TARGET stimulus for information masking experiments
% (Durlach et al., 2003a; Durlach et al., 2003b; Kidd et al., 2002; Kidd et
% al., 2003; Kidd et al., 1998; Oxenham et al., 2003).
%
%

% 2013 Marc van Wanrooij (original code by Peter Bremen)
% e-mail: marcvanwanrooij@neural-code.com

tic

%% Initialization
dispFlag = false;
% if dispFlag
% 	close all
% end

Fs		=	48828;
dBRef	=	94; % dB SPL

if nargin<1
	Par.MaskRand	= true; % Masker timing: true: asynchronous; false: synchronous
	Par.Paradigm	= 2; %  1: random frequencies; 2: constant frequencies across pulses
	Par.ProtectBand	= 1; % Protected Band (in 2013: 0:1/3:2)
	Par.NMaskLow	= 2; % Number of masker tones below target
	Par.NMaskHigh	= 2; % number of masker tones above target
	Par.MaskRng		= 1; % Range (octave) from which maskers can be drawn
Par.NPulse		= 4; % Number of pulses
Par.SigFreq		= 1600; % target frequency


end
%% Parameters
Par.BaseRate	= 10; % Repetition rate of pulses
Par.MaskSPL		= 45; % Masker Level
Par.GaussDur	= 50; % Duration Envelope
Par.GaussSigma	= 6.25; % Sigma Envelope

%% Generation
nPB	=	length(Par.ProtectBand);
for k = 1:nPB
	FreqVec		= pa_oct2bw(200,0:1/6:6+1/6)';
	MaskGainVec	= ones(length(FreqVec),1);
	SigFreq		= Par.SigFreq;
	
	Freq		= createmasker(Par,SigFreq,FreqVec,MaskGainVec);
	
	NMask		= Par.NMaskLow + Par.NMaskHigh;
	myPeriod	= (1000/Par.BaseRate) / NMask;
	if Par.MaskRand
		[~,idx]		=	sort( rand(1,NMask) );
% 		MaskDelayVec		=	myPeriod*(idx-1) + round((rand(1,1)-.5)*10);
				MaskDelayVec		=	myPeriod*(idx-1);

	else
		MaskDelayVec		=	zeros(1,NMask);
	end
	
	FreqMat		= squeeze(Freq);
	AmpMat		= ones(size(FreqMat));
	M			= reconmask(Par,SigFreq,FreqMat,AmpMat,MaskDelayVec,Fs,dispFlag);
	
	SigAmps		= (10^((60-dBRef)/20)) ./ ones(Par.NPulse);
	SigAmps		= ones(size(SigAmps));
	S			= reconsig(Par,SigFreq,SigAmps,0,Fs,dispFlag);
end

M = sum(M,2);
S = pa_appendn(S,length(M)-length(S));


function [FreqMat,AmpMat] = createmasker(Par,SigFreq,FreqVec,MaskGainVec)

% LoadMaskBuffers: For each interval, fill the frequency and amplitude buffers for the info masker

dBRef	=	94;

NMask	=	Par.NMaskLow + Par.NMaskHigh;

Protect	=	SigFreq * 2.^[-Par.ProtectBand Par.ProtectBand];
AmpVec	=	(10^((Par.MaskSPL-dBRef)/20)) ./ MaskGainVec;
AmpVec	=	min(AmpVec, 1.0);								%-- Limit amplitudes to avoid speaker distortion --%

AmpMat	=	zeros(Par.NPulse,NMask);
FreqMat	=	1000 * ones(Par.NPulse,NMask);

for iPulse=1:Par.NPulse
	Fbands		=	getfreqband(FreqVec,Protect,Par.MaskRng);
	
	low			=	find( FreqVec > Fbands(1) & FreqVec < Protect(1) );
	high		=	find( FreqVec < Fbands(2) & FreqVec > Protect(2) );
	
	Nuni		=	0;
	while( Nuni ~= Par.NMaskLow )
		if( length(high) == 1 )
			ldx		=	randi(length(low),Par.NMaskLow+Par.NMaskHigh-1,1);
		else
			ldx		=	randi(length(low),Par.NMaskLow,1);
		end
		Nuni	=	length( unique(ldx) );
	end
	LowFreq		=	FreqVec(low(ldx))';
	LowAmp		=	AmpVec(low(ldx))';
	
	if( length(high) > 1 )
		Nuni		=	0;
		while( Nuni ~= Par.NMaskHigh )
			hdx		=	randi(length(high),Par.NMaskHigh,1);
			Nuni	=	length( unique(hdx) );
		end
	else
		hdx		=	1;
	end
	HighFreq	=	FreqVec(high(hdx))';
	HighAmp		=	AmpVec(high(hdx))';
	
	FreqMat(iPulse,1:NMask)	=	[LowFreq HighFreq];
	AmpMat(iPulse,1:NMask)	=	[LowAmp HighAmp];
end

%-- C: all pulses have the same frequency --%
if( Par.Paradigm == 2 )
	FreqMat		=	repmat(FreqMat(1,:),Par.NPulse,1);
	AmpMat		=	repmat(AmpMat(1,:),Par.NPulse,1);
end

function Fbands = getfreqband(FreqVec,Protect,MaskRng)

Mrng(1)		=	Protect(1) * 2.^-MaskRng;
Mrng(2)		=	Protect(2) * 2.^+MaskRng;

sel			=	FreqVec >= Mrng(1) & FreqVec <= Mrng(2);
Bands		=	FreqVec(sel);
Fbands		=	[Bands(1) Bands(end)];

function M = reconmask(Par,SigFreq,FreqMat,AmpMat,MaskDelayVec,Fs,flag)

if( nargin < 7 )
	flag	=	0;
end

T			=	1 / Par.BaseRate;
Dur			=	Par.GaussDur * 10^-3;
Sigma		=	Par.GaussSigma * 10^-3;
NPulse		=	Par.NPulse;

NPer		=	round( T * Fs );
NSamp		=	round( Dur * Fs ) + 1;
Noct		=	Par.ProtectBand;
ProtectBand	=	[SigFreq / 2^Noct, SigFreq * 2^Noct]  + [-100 100];
CB1_3		=	[SigFreq / 2^(1/6), SigFreq * 2^(1/6)];

NMaskBand	=	Par.NMaskLow + Par.NMaskHigh;
NMaskDelay	=	round( MaskDelayVec*10^-3 * Fs );
NMask		=	NPer*(NPulse+1);

M			=	zeros(NMask,NMaskBand);
start		=	1;
stop		=	start + NSamp - 1;
mtart		=	1;
for k=1:NPulse
	start			=	stop + (NPer - NSamp);
	stop			=	start + NSamp - 1;
	
	for l=1:NMaskBand
		if( NMaskDelay(l) == 0 )
			ttart	=	1;
		else
			ttart	=	NMaskDelay(l);
		end
		ttop		=	ttart + NSamp -1;
		Tmp(ttart:ttop,l)	=	makesine(FreqMat(k,l),AmpMat(k,l),Dur,Sigma,Fs); %#ok<AGROW>
	end
	
	mtop			=	mtart + size(Tmp,1)-1;
	M(mtart:mtop,:)	=	M(mtart:mtop,:) + Tmp;
	mtart			=	start;
end

if( flag )
	M			=	sum(M,2);
	% 	M			=	M ./ max(abs(M));
	
	tm			=	( ( (0:NMask-1) ) * 1/Fs ) * 10^3;
	
	window		=	round( 10^-2*Fs );
	noverlap	=	round( .5*10^-3*Fs );
	Nfft		=	1024;
	
	figure(666)
	subplot(121)
	plot(tm,M,'k-')
	xlim([tm(1) tm(end)])
	ylim([-1.1 1.1])
	xlabel('time [msec]')
	ylabel('amplitude (a.u.)')
	title('masker alone')
	axis('square')
	
	subplot(122)
% 	snd = sum(M,2);
% 	nsamples	= length(snd);
% 	t			= nsamples/Fs*1000;
% 	dt			= 5;
% 	nseg		= t/dt;
% 	segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
% 	noverlap	= round(0.7*segsamples); % 1/3 overlap
% 	window		= segsamples+noverlap; % window size
% 	nfft		= 10000;
% 	spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
% 	% spectrogram(snd)
% 	% colorbar
% 	cax = caxis;
% 	caxis([0.7*cax(1) 1.1*cax(2)])
% 	ylim([100 8000])
% 	set(gca,'YTick',[0.5 1 2 4 8 16]*1000);
% 	set(gca,'YTickLabel',[0.5 1 2 4 8 16]);
% 	% set(gca,'YScale','log');
% 	drawnow
% 	axis square;
% 	caxis([-70 -20])
% 	colormap gray
% 	
% 	% keyboard
		[~,F,T,P] = spectrogram(sum(M,2),window,noverlap,Nfft,Fs);
		T	=	T * 10^3;
		surf(T,F,10*log10(abs(P)+eps),'EdgeColor','none');
		axis tight
		view(0,90);
		hold on
		plot([T(1) T(end)],[ProtectBand(1) ProtectBand(1)],'r-')
		plot([T(1) T(end)],[ProtectBand(2) ProtectBand(2)],'r-')
		plot([T(1) T(end)],[CB1_3(1) CB1_3(1)],'r--')
		plot([T(1) T(end)],[CB1_3(2) CB1_3(2)],'r--')
		set(gca,'YTick',[62.5 125 250 500 1000 2000 4000 8000 16000 32000],'YTickLabel',[62.5 125 250 500 1000 2000 4000 8000 16000 32000])
		set(gca,'YScale','log')
		xlabel('time [ms]')
		ylabel('frequency [Hz]')
		title('masker alone')
		box on
		axis('square')
	caxis([-70 -20]);
		colormap('gray')
end

function S = reconsig(Par,SigFreq,SigAmp,SigDelay,Fs,flag)

if( nargin < 6 )
	flag	=	0;
end

T			=	1 / Par.BaseRate;
Dur			=	Par.GaussDur * 10^-3;
Sigma		=	Par.GaussSigma * 10^-3;
NPulse		=	Par.NPulse;

NPer		=	round( T * Fs );
NSigDelay	=	round( SigDelay*10^-3 * Fs );
NSamp		=	round( Dur * Fs ) + 1;
NSig		=	NPer*NPulse + NSigDelay;
Noct		=	Par.ProtectBand;
ProtectBand	=	[SigFreq / 2^Noct, SigFreq * 2^Noct]  + [-100 100];
CB1_3		=	[SigFreq / 2^(1/6), SigFreq * 2^(1/6)];

S			=	zeros(NSig,1);
start		=	NSigDelay + 1;
stop		=	start + NSamp - 1;
for k=1:NPulse
	S(start:stop,1)	=	makesine(SigFreq,SigAmp(k),Dur,Sigma,Fs);
	start			=	stop + (NPer - NSamp);
	stop			=	start + NSamp - 1;
end

if( flag )
	% 	S			=	S ./ max(abs(S));
	ts			=	( ( (0:NSig-1) ) * 1/Fs ) * 10^3;
	
	window		=	round( 10^-2*Fs );
	noverlap	=	round( .5*10^-3*Fs );
	Nfft		=	1024;
	
	figure(2)
	clf
	
	subplot(2,1,1)
	plot(ts,S,'k-')
	xlim([ts(1) ts(end)])
	ylim([-1.1 1.1])
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('time [msec]')
	ylabel('amplitude [a.u.]')
	title('signal alone')
	axis('square')
	
	subplot(2,1,2)
	[~,F,T,P] = spectrogram(S,window,noverlap,Nfft,Fs);
	T	=	T * 10^3;
	surf(T,F,10*log10(abs(P)+eps),'EdgeColor','none');
	axis tight
	view(0,90);
	hold on
	plot([T(1) T(end)],[ProtectBand(1) ProtectBand(1)],'r-')
	plot([T(1) T(end)],[ProtectBand(2) ProtectBand(2)],'r-')
	plot([T(1) T(end)],[CB1_3(1) CB1_3(1)],'r--')
	plot([T(1) T(end)],[CB1_3(2) CB1_3(2)],'r--')
	set(gca,'YTick',[62.5 125 250 500 1000 2000 4000 8000 16000 32000],'YTickLabel',[62.5 125 250 500 1000 2000 4000 8000 16000 32000])
	set(gca,'YScale','log')
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('time [ms]')
	ylabel('frequency [Hz]')
	title('signal alone')
	box on
	axis('square')
	
	colormap('gray')
end

function [s,x] = makesine(F,G,D,ramp,Fs)

len		=	length(F);

Ndur	=	round( D * Fs );
Nenv	=	round( ramp * Fs );

x		=	0:1/Fs:D;
s		=	nan(len,Ndur+1);
for k=1:len
	s(k,:)	=	G .* sin(2*pi*F(k)*x);
	s(k,:)	=	applyramp(s(k,:),Ndur,Nenv);
end

function out = applyramp(in,N,Ns,flag)

if( nargin < 4 )
	flag	=	0;
end

len		=	length(in)-1;

x		=	(0:N) / len;
Ns		=	Ns / len;
m		=	x(end) / 2;

g		=	exp( -( (x-m) ./ Ns ).^2 );
out		=	g .* in;

if( flag )
	figure(1)
	clf
	
	subplot(2,1,1)
	plot(x,g,'k-')
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('normalized time')
	ylabel('normalized amplitude')
	title('Gaussian filter')
	
	subplot(2,1,2)
	plot(x,out,'k-')
	hold on
	plot(x,in,'r--')
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('normalized time')
	ylabel('normalized amplitude')
	title('original filter (red) & filtered signal (black)')
end
