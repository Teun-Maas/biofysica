function Flow = pa_dtfanalysis(fname)
% PA_DTFANALYSIS(FNAME)
%
% Analyze the acoustic recordings, and produce the directional transfer
% funtion.
%
% See also PA_READHRTF, PA_SWEEP2MAG, PA_READCSV
if nargin<1
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
end

freqrange	= [4000 15000];
locrange	= 90;
NFFT		= 1024;
Nsweep		= 18;
Fs			= 48828.125;
nBegin		= ceil(50/1000*Fs);
D = struct([]);
tic;
for ii = 1:2
	CH = ii;
	
	%% Microphone & Room Acoustics

	cd('E:\MATLAB\PANDA\Audition\HRTF\Knowles\');
	wavfile			= 'snd399.wav';
	wav				= wavread(wavfile);
	


	
	
	%% Limit Location Range

	
	%% Spectrum
	cd(['E:\DATA\Head Related Transfer Functions\' fname(1:5) '\' fname(1:16)]);
	hrtf			= pa_readhrtf(fname);
	
	csvfile			= pa_fcheckext(fname,'csv');
	[exp,chan,mlog] = pa_readcsv(csvfile);
	sel             = mlog(:,5) == 2; %SND1 is the Schroeder Sweep, SND2 is the zero-Sweep
	Theta           = mlog(sel,6);
	Phi             = mlog(sel,7);
	[Azimuth,Elevation] = pa_fart2azel(Theta,Phi);
	Azimuth         = round(Azimuth*2)/2; % We can round to integers
	Elevation       = round(Elevation*2)/2;
	[Spec,freq,Azimuth,Elevation] = pa_sweep2mag(hrtf,wav,Azimuth,Elevation,CH,NFFT,Nsweep,Fs,nBegin,0);
	toc
	
	
		sel = Elevation==0;
	A				= Azimuth(sel);
	E				= Elevation(sel);
	
	azSpec			= Spec(:,sel);
	Mag				= azSpec;
	
	%% Remove Directional Component
	rms                 =   sqrt(mean(Mag.^2,2));
	rms                 =   repmat(rms,1,size(Mag,2));
	DTF                 =   Mag./rms;
	DTF(1,:)			= 0;
	
	%% Remove echo/reflections (low-frequency)
	DTF                 =   pa_remecho(DTF,Fs);
	
	%% Smooth
	DTF					= pa_locsmooth(DTF,2,A,E);
	DTF                 = pa_hsmooth(DTF,6);
		
	
	D(ii).dtf = DTF;
end

DTF1 = D(1).dtf;
DTF2 = D(2).dtf;

keyboard
%%
	col				= pa_statcolor(ceil(size(DTF2,2)/2)*2,[],[],[],'def',1);

Az = -A;
figure(1)
set(gcf,'DefaultAxesColorOrder',col)
% subplot(221)
% 
% h = plot(freq(1:513),DTF2)
% box off
% x = pa_oct2bw(500,-1:4);
% set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log');
% xlim(pa_oct2bw(500,[2 5]))
% axis square;
% colorbar
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% title('DTF');
% % legend(h,num2str(Az));
% pa_text(0.01,0.9,char(64+1));

subplot(131)
colimg				= pa_statcolor(ceil(size(DTF2,2)/2)*2,[],[],[],'def',8);
colormap(colimg)
contourf(freq(1:513),Az,DTF2',-100:20)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5)
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[2 5]))
ylabel('Azimuth (deg)');
xlabel('Frequency (Hz)');
title('DTF');
pa_text(0.01,0.9,char(64+2));

sel = freq(1:513)>2000 & freq(1:513)<12000;
muHSE = mean(DTF2(sel,:));
subplot(132)
plot(Az,muHSE,'k-');
axis square;
box off
colorbar
xlabel('Azimuth (deg)');
ylabel('Power (dB)');
title({'Average pinna/HSE';'2 kHz < Freq < 12 kHz'});
ylim([-15 15])
pa_horline(0,'k:');
pa_verline([-30 30],'k:');
pa_text(0.01,0.9,char(64+3));

subplot(133)
DTF3 = bsxfun(@minus,DTF2,muHSE);
colimg				= pa_statcolor(ceil(size(DTF3,2)/2)*2,[],[],[],'def',8);
colormap(colimg)
contourf(freq(1:513),-Az,flipud(DTF3'),-20:20)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[2 5]))
xlabel('Frequency (Hz)');
ylabel('Azimuth (deg)');
title('Frequency-dependent pinna/HSE');
pa_text(0.01,0.9,char(64+4));

% return
%%
close all
figure(2)
col				= pa_statcolor(513,[],[],[],'def',1);
a = Az;

for ii = 50:10:400
	indx = ii;
	sel = freq==freq(indx);
	dtf = DTF2(sel,:);
	whos freq DTF
	subplot(121)
	plot(a,dtf,'k-','Color',col(indx,:));
	hold on
end
box off
axis square;
set(gca,'TickDir','out');
ylim([-40 20]);
ylabel('HSE (dB)');
xlabel('Azimuth (deg)');
title(round(freq(indx)));

% %% HSE = 0
for ii = 150
	indx = ii;
	sel = freq==freq(indx);
	dtf = DTF2(sel,:); % function of alpha
	plot(a,dtf,'k-','Color',col(indx,:));
	naz = numel(a);
	udtf = -30:30;
	ndtf = numel(udtf);
	r = NaN(naz,ndtf);
	for azIdx = 1:naz
	r(azIdx,:) = normpdf(udtf,+dtf(azIdx),11);
	end
% 	whos r a
% 	plot(a,mean(r,2))
subplot(132)
	imagesc(Az,udtf,r')
	hold on

	subplot(133)
	sel = udtf == 0;
	plot(Az,r(:,sel),'k-')
	hold on

end
for ii = 1:3
	subplot(1,3,ii)
box off
axis square;
set(gca,'TickDir','out');
set(gca,'YDir','normal');
end
% ylim([-40 20]);
% ylabel('HSE (dB)');
% xlabel('Azimuth (deg)');
% title(round(freq(indx)));


