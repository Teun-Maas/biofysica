function pa_headshadoweffect_lv(fname)
% PA_DTFANALYSIS(FNAME)
%
% Analyze the acoustic recordings, and produce the directional transfer
% funtion.
%
% See also PA_READHRTF, PA_SWEEP2MAG, PA_READCSV
if nargin<1
	close all
% 	clear all
% 		fname				= 'SB-PB-2008-02-20-0004.hrtf';
	% 	fname				= 'MW-BG-2010-05-06-0001.hrtf';
	% 	fname				= 'MW-MA-2010-07-02-0001.hrtf';
% 		fname				= 'MW-RW-2010-07-02-0002.hrtf';
	% 	fname				= 'MW-BG-2010-07-14-0002.hrtf';
	% 	fname				= 'MW-BS-2010-07-14-0004.hrtf';
% 		fname				= 'MW-XX-2010-07-15-0004.hrtf';
% 		fname				= 'MA-DE-2010-07-19-0001.hrtf';
% 		fname				= 'MA-OT-2010-08-04-0004.hrtf';
	fname				= 'DM-TH-2012-03-14-0001.hrtf';
end

NFFT		= 1024;
Nsweep		= 18;
Fs			= 48828.125;
nBegin		= ceil(50/1000*Fs);

tic;
D(2).dtf = 0;
for ii = 1:2
	CH = ii;
	
	%% Microphone & Room Acoustics	
	cd('/Users/marcw/matlab/PandA/Audition/HRTF/Knowles');
	wavfile			= 'snd399.wav';
	wav				= audioread(wavfile);
	
	%% Spectrum
	cd(['/Users/marcw/DATA/Head Related Transfer Functions/' fname(1:5) '/' fname(1:16)]);
	hrtf			= pa_readhrtf(fname);
	csvfile			= pa_fcheckext(fname,'csv');
	[~,~,mlog] = pa_readcsv(csvfile);
	sel             = mlog(:,5) == 2; %SND1 is the Schroeder Sweep, SND2 is the zero-Sweep
	Theta           = mlog(sel,6);
	Phi             = mlog(sel,7);
	[Azimuth,Elevation] = pa_fart2azel(Theta,Phi);
	Azimuth         = round(Azimuth*2)/2; % We can round to integers
	Elevation       = round(Elevation*2)/2;
	[Spec,freq,Azimuth,Elevation] = pa_sweep2mag(hrtf,wav,Azimuth,Elevation,CH,NFFT,Nsweep,Fs,nBegin,0);

	sel				= Elevation==0;
	A				= Azimuth(sel);
	E				= Elevation(sel);
	azSpec = Spec(:,sel);
	Mag				= azSpec;
	%% Remove Directional Component
	rms                 =   sqrt(mean(Mag.^2,2));
	DTF                 =   bsxfun(@rdivide,Mag,rms);
	DTF(1,:)			= 0;
	
	%% Remove echo/reflections (low-frequency)
	DTF                 =   pa_remecho(DTF,Fs);
	
	%% Smooth
	DTF					= pa_locsmooth(DTF,10,A,E);
	DTF                 = pa_hsmooth(DTF,6);
	col					= pa_statcolor(size(DTF,2),[],[],[],'def',1);
	
	
	colimg				= pa_statcolor(ceil(size(DTF,2)/2)*2,[],[],[],'def',8);
	colormap(colimg)
	
	%%
	figure(1)
	set(gcf,'DefaultAxesColorOrder',col)
	subplot(121)
	plot(freq(1:513),DTF)
	box off
	x = pa_oct2bw(500,-1:4);
	set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log');
	xlim(pa_oct2bw(500,[0 5]))
	axis square;
	
	figure(3)
	subplot(132)
	contourf(freq(1:513),-A,DTF',-100:20)
	shading flat
	box off
	axis square
	set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log','YDir','normal');
	colorbar
	caxis([-15 5])
	xlim(pa_oct2bw(500,[0 5]))
	
	
	
	D(ii).dtf = DTF;
end

xlabel('Frequency (Hz)');
ylabel('Azimuth (deg)');
title({'Overview';'Power (dB(V))'});
DTF1 = D(1).dtf;
DTF2 = D(2).dtf;

%%
figure;
% warning off
Az	= -A;
figure(3)
set(gcf,'DefaultAxesColorOrder',col);
subplot(131);

sel = abs(Az)<90;
X	= Az(sel);
Y	= DTF2(:,sel)';
Y	= bsxfun(@minus,Y,Y(19,:));
t = Y(:,50:10:end);
% t = bsxfun(@plus,t,1:size(t,2))-1;
t = t';
plot(X,t);
hold on
drawnow

B		= NaN(512,1);
Xnew	= -90:90;
Ynew	= NaN(numel(freq),numel(Xnew));
beta0	= [9.7 0.02 0.27 0];
fun		= @(beta,x)(beta(1)*sin(beta(2)*x+beta(3))+beta(4));
for ii	= 1:512
	idx				= ii-10:ii+10;
	idx(idx<1)		= 1;
	idx(idx>512)	= 512;
	idx				= unique(idx);
	y				= nanmean(Y(:,idx),2);
	beta			= nlinfit(X,y,fun,beta0);
	mdl				= fitnlm(X,y,fun,beta0);
	[ynew,ynewci]	= predict(mdl,Xnew');
	Ynew(ii,:)		= ynew;
	B(ii)			= beta(1);
end


t = Ynew(50:30:end,:);
t = bsxfun(@plus,t',1:size(t,1))-1;
t = t';
subplot(131)
% plot(Xnew,t,'Color',[.7 .7 .7]);
box off
ylim([-15 17])
xlim([-100 100])
pa_horline(0,'k-');
set(gca,'XTick',-90:30:90,'TickDir','out');
box off
axis square;
xlabel('Azimuth (deg)');
ylabel('Power (dB(V))');
title({'Head shadow and pinna reflections';'for various frequencies'})

subplot(133)
semilogx(freq(1:512),B,'k-')
ylim([0 12])
xlim([500 20000])
set(gca,'XTick',pa_oct2bw(1000,0:4),'XTickLabel',pa_oct2bw(1000,0:4)/1000,'TickDir','out');
box off
axis square;
xlabel('Frequency (kHz)');
ylabel('Fitted sine amplitude (dB(V))');
title('Fitted head shadow as a function of frequency');
drawnow

pa_datadir;
print('-depsc','-painters','head shadow');


%%
