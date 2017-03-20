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

freqrange	= [2000 15000];
locrange	= 90;
NFFT		= 1024;
Nsweep		= 18;
Fs			= 48828.125;
nBegin		= ceil(50/1000*Fs);

tic;
for ii = 1:2
	CH = ii;
	
	%% Microphone & Room Acoustics

	cd('/Users/marcw/Gitlab/panda/Audition/HRTF/Knowles');
	wavfile			= 'snd399.wav';
	wav				= audioread(wavfile);
	


	
	
	%% Limit Location Range

	
	%% Spectrum
	cd(['/Users/marcw/DATA/Head Related Transfer Functions/' fname(1:5) '/' fname(1:16)]);
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
	A				= Azimuth;
	E				= Elevation;
	
	azSpec = Spec;
	Mag				= azSpec;
	whos A E Mag
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
	col					= pa_statcolor(size(DTF,2),[],[],[],'def',1);
	
	
	colimg				= pa_statcolor(ceil(size(DTF,2)/2)*2,[],[],[],'def',8);
	colormap(colimg)
	
	%%
	figure(1)
	set(gcf,'DefaultAxesColorOrder',col)
	subplot(121)
	plot(freq(1:513),DTF)
	box off
	x = pa_oct2bw(500,-1:4)
	set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log');
	xlim(pa_oct2bw(500,[0 5]))
	axis square;
	
	subplot(122)
	contourf(freq(1:513),A,DTF',-100:20)
	shading flat
	box off
	axis square
	set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log','YDir','normal');
	colorbar
	caxis([-15 5])
	xlim(pa_oct2bw(500,[0 5]))
	
	%%
	
	
	
	D(ii).dtf = DTF;
end

DTF1 = D(1).dtf;
DTF2 = D(2).dtf;


%%
Az = -A;
figure(3)
set(gcf,'DefaultAxesColorOrder',col)
subplot(221)

h = plot(freq(1:513),DTF2)
box off
x = pa_oct2bw(500,-1:4);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log');
xlim(pa_oct2bw(500,[2 5]))
axis square;
colorbar
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('DTF');
% legend(h,num2str(Az));
pa_text(0.01,0.9,char(64+1));

subplot(222)
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
subplot(223)
plot(Az,muHSE,'ko-','MarkerFaceColor','w');
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

subplot(224)
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

% set(gcf,'DefaultAxesColorOrder',col)
% 
% h = plot(freq(1:513),DTF3)
% box off
% x = pa_oct2bw(500,-1:4);
% set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','log');
% xlim(pa_oct2bw(500,[2 5]))
% axis square;
% colorbar
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% title('DTF');
% legend(h,num2str(A));
pa_datadir
print('-depsc','-painter',mfilename);
%%
% keyboard
%%
DTF = DTF2-DTF1;
DTF = DTF2;

el = E;
az = A;
sel = el==0;
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
figure(4)
set(gcf,'DefaultAxesColorOrder',col)
subplot(241)
colimg				= pa_statcolor(ceil(size(dtf,2)/2)*2,[],[],[],'def',8);
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 5]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(245)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar

el = E;
az = A;
[az,el] = pa_2drotate(az,el,-30);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(4)
set(gcf,'DefaultAxesColorOrder',col)
subplot(242)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 5]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(246)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar


el = E;
az = A;
[az,el] = pa_2drotate(az,el,-60);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(4)
set(gcf,'DefaultAxesColorOrder',col)
subplot(243)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 5]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(247)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar

el = E;
az = A;
[az,el] = pa_2drotate(az,el,-90);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(4)
set(gcf,'DefaultAxesColorOrder',col)
subplot(244)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 5]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(248)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar

pa_datadir
print('-depsc','-painter',[mfilename '1']);

%%

%%
% close all
warning off
DTF = DTF2-DTF1;

el = E;
az = A;
sel = el==0;
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
figure(5)
set(gcf,'DefaultAxesColorOrder',col)
subplot(241)
colimg				= pa_statcolor(ceil(size(dtf,2)/2)*2,[],[],[],'def',8);
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:3);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 5]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(245)
sel = az==0;
HH = dtf(:,sel);
plot(freq(1:513),dtf(:,sel),'k-')
box off
axis square
x = pa_oct2bw(500,0:5);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear');
xlim([2000 16000]);
% axis([-90 90 -20 20]);
colorbar
R = NaN(size(dtf,2),1);
for ii = 1:size(dtf,2);
	r = corrcoef(HH,dtf(:,ii));
	R(ii) = r(2);
end
figure(6)
subplot(141);
plot(az,R,'k-');
axis square
box off
wa = sum(az.*R)/sum(R);
pa_verline(wa)
title(round(wa))


el = E;
az = A;
[az,el] = pa_2drotate(az,el,-30);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(5)
set(gcf,'DefaultAxesColorOrder',col)
subplot(242)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:3);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 3]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(246)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar
R = NaN(size(dtf,2),1);
for ii = 1:size(dtf,2);
	r = corrcoef(HH,dtf(:,ii));
	R(ii) = r(2);
end
figure(6)
subplot(142);
plot(az,R,'k-');
axis square
box off
wa = sum(az.*R)/sum(R);
pa_verline(wa)
title(round(wa))

el = E;
az = A;
[az,el] = pa_2drotate(az,el,-60);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(5)
set(gcf,'DefaultAxesColorOrder',col)
subplot(243)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:3);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 3]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(247)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar
R = NaN(size(dtf,2),1);
for ii = 1:size(dtf,2);
	r = corrcoef(HH,dtf(:,ii));
	R(ii) = r(2);
end
figure(6)
subplot(143);
plot(az,R,'k-');
axis square
box off
wa = sum(az.*R)/sum(R);
pa_verline(wa)
title(round(wa))

el = E;
az = A;
[az,el] = pa_2drotate(az,el,-90);
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
sum(sel)
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
% close all
figure(5)
set(gcf,'DefaultAxesColorOrder',col)
subplot(244)
colormap(colimg)
contourf(freq(1:513),az,dtf',-50:50)
shading flat
box off
axis square
x = pa_oct2bw(500,0:3);
set(gca,'TickDir','out','XTick',x,'XtickLabel',x,'XScale','linear','YDir','normal');
colorbar
caxis([-10 10])
xlim(pa_oct2bw(500,[0 3]))
title('ILD');
ylabel('Azimuth (deg)');
subplot(248)
plot(az,mean(dtf),'k-')
box off
axis square
axis([-90 90 -20 20]);
colorbar
R = NaN(size(dtf,2),1);
for ii = 1:size(dtf,2);
	r = corrcoef(HH,dtf(:,ii));
	R(ii) = r(2);
end
figure(6)
subplot(144);
plot(az,R,'k-');
axis square
box off
wa = sum(az.*R)/sum(R);
pa_verline(wa)
title(round(wa))


%%
DTF = DTF2-DTF1;

el = E;
az = A;
sel = el==0;
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);
sel = az==0;
HH = dtf(:,sel);


% close all
roll = -120:5:120;
W = NaN(size(roll));
for jj = 1:numel(roll);
	
el = E;
az = A;
[az,el] = pa_2drotate(az,el,-1*roll(jj));
az = round(az/2)*2;
el = round(el/2)*2;
sel = el==0;
el = el(sel);
az = az(sel);
dtf = DTF(:,sel);

R = NaN(size(dtf,2),1);
for ii = 1:size(dtf,2);
	r = corrcoef(HH,dtf(:,ii));
	R(ii) = r(2);
end
wa = nansum(az.*R)/nansum(R);
W(jj) = wa;
end

figure(7)
indx = ismember(roll,-120:30:120);
W = smooth(W,10);
plot(roll(indx),W(indx)*15/37,'k-','LineWidth',2);
hold on
% plot(roll,W,'ko','MarkerFaceColor','w');

axis square
box off
xlim([-130 130]);
ylim([-20 20]);
set(gca,'TickDir','out');
xlabel('Roll angle (deg)');
ylabel('Level difference (dB)');
pa_horline(0,'k:');
pa_verline(0,'k:');

pa_datadir
print('-depsc','-painter',[mfilename '2']);


%%

% keyboard

%
% %% Left and Right Ear Averaged
% DTF = (D(1).dtf+D(2).dtf)/2;
% figure
% subplot(121)
% sel = E==-40;
% self = F>3500;
% [mn1,indx1] = min(DTF(self,sel));
% semilogx(F,DTF(:,sel),'r-');
% hold on
% f = F(self);
% plot(f(indx1),mn1,'ko');
% Flow = f(indx1);
%
% sel = E==40;
% [mn2,indx2] = min(DTF(:,sel));
% semilogx(F,DTF(:,sel),'b-');
% plot(F(indx2),mn2,'ko');
% x = [1000 2000 4000 8000 16000];
% set(gca,'Xtick',x,'XTickLabel',x/1000);
% xlim(minmax(F))
% axis square
% xlabel('Frequency (kHz)');
% ylabel('Power (dB)');
% title(f(indx1));
%
% maxdtf  = max(max(DTF));
% mindtf  = min(min(DTF));
% cnt = mindtf:1:maxdtf;
% subplot(122)
% contourf(F,uE,DTF',20);
% hold on
% plot(f(indx1),-40,'wo','LineWidth',3);
% plot(F(indx2),40,'wo','LineWidth',3);
% x = [1000 2000 4000 8000 16000];
% set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
% xlim(minmax(F))
% axis square
% xlabel('Frequency (kHz)');
% ylabel('Elevation (deg)');
% title(F(indx2));
%
% figure
% maxdtf  = max(max(DTF));
% mindtf  = min(min(DTF));
% cnt = mindtf:1:maxdtf;
% contourf(F,uE,DTF',20);
% hold on
% x = [1000 2000 4000 8000 16000];
% set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
% xlim(minmax(F))
% axis square
% xlabel('Frequentie (kHz)');
% ylabel('Vertikale richting (graden)');
