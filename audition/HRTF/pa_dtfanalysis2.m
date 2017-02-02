% function Flow = pa_dtfanalysis(fname)
% PA_DTFANALYSIS(FNAME)
%
% Analyze the acoustic recordings, and produce the directional transfer
% funtion.
%
% See also PA_READHRTF, PA_SWEEP2MAG, PA_READCSV
% if nargin<1
	close all
	clear all
% 	fname				= 'SB-PB-2008-02-20-0004.hrtf';
% 	fname				= 'MW-BG-2010-05-06-0001.hrtf';
% 	fname				= 'MW-MA-2010-07-02-0001.hrtf';
% 	fname				= 'MW-RW-2010-07-02-0002.hrtf';
% 	fname				= 'MW-BG-2010-07-14-0002.hrtf';
% 	fname				= 'MW-BS-2010-07-14-0004.hrtf';
% 	fname				= 'MW-XX-2010-07-15-0004.hrtf';
% 	fname				= 'MA-DE-2010-07-19-0001.hrtf';
% 	fname				= 'MA-OT-2010-08-04-0004.hrtf';
	fname				= 'DM-TH-2012-03-14-0001.hrtf';
	
	
% end

freqrange	= [2000 15000];
locrange	= 90;
NFFT		= 1024;
Nsweep		= 18;
Fs			= 48828.125;
nBegin		= ceil(50/1000*Fs);

tic;
for ii = 1:2
	CH = ii;
	
	cd('E:\MATLAB\PANDA\Audition\HRTF\Knowles\');
	wavfile			= 'snd399.wav';
	wav				= wavread(wavfile);
	%% Spectrum
	cd(['E:\DATA\Head Related Transfer Functions\' fname(1:5) '\' fname(1:16)]);
	hrtf			= pa_readhrtf(fname);
	csvfile			= pa_fcheckext(fname,'csv');
	[exp,chan,mlog] = pa_readcsv(csvfile);
	sel             = mlog(:,5) == 2; %SND1 is the Schroeder Sweep, SND2 is the zero-Sweep
	Theta           = mlog(sel,6);
	Phi             = mlog(sel,7);
	[Azimuth,Elevation] = pa_fart2azel(Theta,Phi);
% 	Azimuth         = round(Azimuth*2)/2; % We can round to integers
% 	Elevation       = round(Elevation*2)/2;
% 	plot(Azimuth,Elevation,'ko','MarkerFaceColor','w')
	
	%% Limit Location Range
	sel				= Azimuth==0 & abs(Elevation)<locrange;
	Azimuth			= Azimuth(sel);
	Elevation		= Elevation(sel);
	hrtf			= hrtf(:,:,sel);
	[Spec,freq,Azimuth,Elevation] = pa_sweep2mag(hrtf,wav,Azimuth,Elevation,CH,NFFT,Nsweep,Fs,nBegin,0);
	toc
	

	
	Mag				= Spec;
	
	%% Remove Directional Component
	rms                 =   sqrt(mean(Mag.^2,2));
	rms                 =   repmat(rms,1,size(Mag,2));
	DTF                 =   Mag./rms;
	DTF(1,:)			= 0;
	
	%% Remove echo/reflections (low-frequency)
	DTF                 =   pa_remecho(DTF,Fs);
	toc
	
	%% Smooth
	DTF					= pa_locsmooth(DTF,2,Azimuth,Elevation);
	DTF                 = pa_hsmooth(DTF,6);
	toc
	
	%% Select Frequency
	F                   =   Fs/2*(0:size(DTF,1)-1)/size(DTF,1);
	sel                 =   F>freqrange(1) & F<freqrange(2);
	F                   =   F(sel);
	DTF                 =   DTF(sel,:);
	
	uE					= unique(Elevation);
	
	tmp					= NaN(size(DTF,1),length(uE));
	for jj				= 1:length(uE)
		sel				= Elevation == uE(jj);
		mu				= mean(DTF(:,sel),2);
		tmp(:,jj)		= mu;
	end
	DTF = tmp;
	
	%% Graph
	figure
	subplot(121)
	for jj = 1:size(DTF,2)
		h1 = semilogx(F,DTF(:,jj)+2*jj,'k-');
		set(h1,'Color',[.8 .8 .8]);
		hold on
	end
	% set(h2,'LineWidth',3);
	x = [1000 2000 4000 8000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000);
	xlim(minmax(F))
	axis square
xlabel('Frequency (kHz)');
ylabel('Power (dB) + elevation offset');
	
	maxdtf  = max(max(DTF));
	mindtf  = min(min(DTF));
	cnt = linspace(mindtf,maxdtf,50);
	
	subplot(122)
	contourf(F,uE,DTF',cnt);
	x = [1000 2000 4000 8000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	xlim(minmax(F))
	axis square
	shading flat
xlabel('Frequency (kHz)');
ylabel('Elevation (deg)');

	D(ii).dtf = DTF;
end

%% Left and Right Ear Averaged
DTF = (D(1).dtf+D(2).dtf)/2;
figure
subplot(121)
sel = Elevation==-40;
self = F>3500;
[mn1,indx1] = min(DTF(self,sel));
semilogx(F,DTF(:,sel),'r-');
hold on
f = F(self);
plot(f(indx1),mn1,'ko');
Flow = f(indx1);

sel = Elevation==40;
[mn2,indx2] = min(DTF(:,sel));
semilogx(F,DTF(:,sel),'b-');
plot(F(indx2),mn2,'ko');
x = [1000 2000 4000 8000 16000];
set(gca,'Xtick',x,'XTickLabel',x/1000);
xlim(minmax(F))
axis square
xlabel('Frequency (kHz)');
ylabel('Power (dB)');
title(f(indx1));

maxdtf  = max(max(DTF));
mindtf  = min(min(DTF));
cnt = mindtf:1:maxdtf;
subplot(122)
contourf(F,uE,DTF',20);
hold on
plot(f(indx1),-40,'wo','LineWidth',3);
plot(F(indx2),40,'wo','LineWidth',3);
x = [1000 2000 4000 8000 16000];
set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
xlim(minmax(F))
axis square
xlabel('Frequency (kHz)');
ylabel('Elevation (deg)');
title(F(indx2));

figure
maxdtf  = max(max(DTF));
mindtf  = min(min(DTF));
cnt = mindtf:1:maxdtf;
contourf(F,uE,DTF',20);
hold on
x = [1000 2000 4000 8000 16000];
set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
xlim(minmax(F))
axis square
xlabel('Frequentie (kHz)');
ylabel('Vertikale richting (graden)');
