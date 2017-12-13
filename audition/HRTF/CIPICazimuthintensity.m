clear all
close all
clc

%% Our Data
% Flow1 = pa_dtfanalysis('MW-BS-2010-07-14-0004.hrtf'); age1 = 65;
% Flow2 = dtfanalysis('MW-XX-2010-07-15-0004.hrtf'); age2 = 80;
% Flow3 = dtfanalysis('MW-X2-2010-07-15-0004.hrtf'); age3 = 76;
% Flow4 = dtfanalysis('MA-DE-2010-07-19-0001.hrtf'); age4 = 65;
% Flow5 = dtfanalysis('MA-OT-2010-08-04-0004.hrtf'); age5 = 73;

%% Anthropometry
c = '/Users/marcw/Gitlab/thirdparty/CIPIC/anthropometry';
cd(c)
load anthro
whos
M = zeros(5,24);
HS = M;
nid = length(id);
for ii = 1:nid
	sub = id(ii);
	sel		= id==sub;
	earsize1 = sum(D(sel,[1 2]),2);
	earsize2 = sum(D(sel,[9 10]),2);
	% 	earsize = D(sel,5);
	%%  Initialize globals
	start_flag = 0;             % Prevents any action until data is available
	log_scale  = 0;             % Default to use a logarithmic frequency scale
	do_smooth  = 1;             % Default to use constant-Q smoothing
	Q          = 8;             % For constant-Q smoothing, of course
	fs         = 44100;         % Sampling frequency in Hz
	fNy        = fs/2000;       % Nyquist frequency in kHz
	Ia         = 13;            % Initial azimuth index   (0 deg)
	Ie         = 9;             % Initial elevation index (0 deg)
	
	A = [-80 -65 -55 -45:5:45 55 65 80];   LA = length(A)  % Azimuths
	E = -45:(360/64):235;                  LE = length(E);  % Elevation
	LT = 200;                                               % Number of time points
	T = (0:(1/fs):((LT-1)/fs))*1000;                        % Times in ms
	subject_numbers = sub_nums; % Load cell array of subject numbers from file
	NFFT = LT;                  % More freq resolution if NFFT > LT, but slows computing
	
	Fmin = 0.5; Fmax = 16;      % Min and max frequency in kHz
	Tmin = 0.3; Tmax = 2.5;     % Min and max time in ms
	miI =  0.0; mxI = 0.75;     % Min and max ITD in ms
	miT = max(find(T<=Tmin));   % Min time index
	mxT = min(find(T>=Tmax));   % Max time index
	dBmin = -30;                % Lower limit for spectral plots in dB
	dBmax = +20;                % Upper limit for spectral plots in dB
	
	%% Data
	
	if sub<10
		c = ['/Users/marcw/Gitlab/thirdparty/CIPIC/standard_hrir_database/subject_00' num2str(sub)];
	elseif sub<100
		c = ['/Users/marcw/Gitlab/thirdparty/CIPIC/standard_hrir_database/subject_0' num2str(sub)];
	elseif sub<1000
		c = ['/Users/marcw/Gitlab/thirdparty/CIPIC/standard_hrir_database/subject_' num2str(sub)];
	end
	cd(c);
	load hrir_final
	hl		= squeeze(hrir_l(:,Ie,:))';               % Normalize data [-1 1]
	hr		= squeeze(hrir_r(:,Ie,:))';
	
	[AL]	= freq_resp(hl,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	[AR,F] = freq_resp(hr,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	LF		= length(F);
	FINIT	= [AL dBmax*ones(LF,1) AR];
	
	sel		= abs(E)<90;
	E		= E(sel);
	A		= A(sel);
	AL		= AL(:,sel);
	AR		= AR(:,sel);
	
	
	figure(1)
	clf
	subplot(131);
	contourf(F,A,AL',20)
	hold on
	x = [500 1000 2000 4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	nicegraph
xlabel('Frequency (kHz)');
ylabel('Head Shadow Left (dB)');

	
	subplot(132);
	contourf(F,A,AR',20);
	hold on;
	x = [500 1000 2000 4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	nicegraph
ylabel('Head Shadow Right (dB)');

	ILD = AR'-AL';
% 	ILD = (AR'+fliplr(AL'))/2;
% 	ILD = AR';
xlabel('Frequency (kHz)');

	subplot(133);
	contourf(F,A,ILD,20);
	hold on;
	x = [500 1000 2000 4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	nicegraph
xlabel('Frequency (kHz)');
ylabel('ILD (dB)');

	%%
	f = [500 1000 2000 4000 8000];
	nf = numel(f);
	for fIdx = 1:nf
	[~,idx]	= min(abs(F-f(fIdx)));
	m = ILD(:,idx)'/nid;
	hs = AR(idx,:)/nid;
% 	m = smooth(m',10)';
	whos m M
	M(fIdx,:) = m+M(fIdx,:);
		HS(fIdx,:) = hs+HS(fIdx,:);

	end
	figure(2)
	plot(A,ILD(:,idx),'-');
	hold on
	
% 	drawnow
	
end
% SPKbubblegraph(x,y);
Mori = M;
%%


%%
M = Mori;
sel = abs(A)<90;
figure(3)
clf
subplot(121)
plot(A(sel),HS(:,sel),'o-','LineWidth',2,'MarkerFaceColor','w');
axis square
box off
ylabel('Head shadow (dB)');
xlabel('Azimuth (deg)');
xlim([-95 95]);
ylim([-30 30]);
set(gca,'TickDir','out',...
	'XTick',-90:30:90,'YTick',-25:5:25);
horline(0,'k:');
h = text([0 0 0 0 0]+65,HS(:,end), num2str(f'/1000));
set(h,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',12)
h = text([0 0 0 0 0]-85,HS(:,1), num2str(f'/1000));
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
title('CIPIC Database');
set(gca,'FontSize',15)
nicegraph

subplot(122)
plot(A(sel),M(:,sel),'o-','LineWidth',2,'MarkerFaceColor','w');
axis square
box off
ylabel('ILD (dB)');
xlabel('Azimuth (deg)');
xlim([-95 95]);
ylim([-30 30]);
set(gca,'TickDir','out',...
	'XTick',-90:30:90,'YTick',-25:5:25);
horline(0,'k:');
h = text([0 0 0 0 0]+65,M(:,end), num2str(f'/1000));
set(h,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',12)
h = text([0 0 0 0 0]-85,M(:,1), num2str(f'/1000));
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
title('CIPIC Database');
set(gca,'FontSize',15)
nicegraph

figure(1)
% savegraph('HS','png');
savegraph('HS','eps');
% savegraph('HS','pdf');
