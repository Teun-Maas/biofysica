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
c = '/Users/marcw/MATLAB/Third Party/CIPIC/anthropometry';
cd(c)
load anthro
whos
M = [];
for ii = 1:length(id)
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
	
	A = [-80 -65 -55 -45:5:45 55 65 80];   LA = length(A);  % Azimuths
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
		c = ['/Users/marcw/MATLAB/Third Party/CIPIC/standard_hrir_database/subject_00' num2str(sub)];
	elseif sub<100
		c = ['/Users/marcw/MATLAB/Third Party/CIPIC/standard_hrir_database/subject_0' num2str(sub)];
	elseif sub<1000
		c = ['/Users/marcw/MATLAB/Third Party/CIPIC/standard_hrir_database/subject_' num2str(sub)];
	end
	cd(c);
	load hrir_final
	hl = squeeze(hrir_l(Ia,:,:))';               % Normalize data [-1 1]
	hr = squeeze(hrir_r(Ia,:,:))';
	
	[AL,Fr] = freq_resp(hl,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	[AR,Fr] = freq_resp(hr,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	LF		= length(Fr);
	FINIT = [AL dBmax*ones(LF,1) AR];
	
	sel = abs(E)<90;
	E = E(sel);
	AL = AL(:,sel);
	AR = AR(:,sel);
	sel = Fr>=3800 & Fr<=15000;
	F = Fr(sel);
	AL = AL(sel,:);
	AR = AR(sel,:);
	sel = F<7500;
	[mn1,indx1] = nanmin(AL(sel,2));
	
	sel = F<7500;
	[mx2,indx2] = nanmin(AR(sel,2));
	
	sel = F<14000;
	[mn3,indx3] = nanmin(AL(sel,16));
	
	sel = F<14000;
	[mx4,indx4] = nanmin(AR(sel,16));
	
	
	figure(1)
	clf
	subplot(221);
	contourf(F,E,AL',20)
	hold on
	x = [4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	title(earsize1)
	axis square;
	plot(F(indx3),E(16),'wo','LineWidth',2);
	plot(F(indx1),E(2),'wo','LineWidth',2);
	
	subplot(223)
	plot(F,AL(:,16),'ko-','MarkerFaceColor','w')
	hold on
	% 	plot(F,AL(:,16),'bo-','MarkerFaceColor','w')
	plot(F,AL(:,2),'ro-','MarkerFaceColor','w')
	verline(F(indx3));
	verline(F(indx1),'r--');
	x = [4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	xlim([min(x) max(x)]);
	horline;
	
	subplot(222);
	contourf(F,E,AR',20);
	
	hold on;
	x = [4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	title(F(indx1));
	plot(F(indx4),E(16),'wo','LineWidth',2);
	plot(F(indx2),E(2),'wo','LineWidth',2);
	
	subplot(224)
	plot(F,gradient(AR(:,16)),'ko-','MarkerFaceColor','w')
	hold on
	plot(F,gradient(AR(:,2)),'ro-','MarkerFaceColor','w')
	verline(F(indx4));
	verline(F(indx2),'r--');
	x = [4000 6000 8000 10000 12000 16000];
	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	axis square;
	xlim([min(x) max(x)]);
	horline;
	
	% 	subplot(212);
	% 	contourf(F,E,AR'+AL',20)
	% 	x = [4000 6000 8000 10000 12000 16000];
	% 	set(gca,'Xtick',x,'XTickLabel',x/1000,'XScale','log');
	% 	axis square;
	
	drawnow
	% 	return
	pause
	
	M = [M;[earsize1 earsize2 F(indx1) F(indx2) F(indx3) F(indx4)]];
end

%%
figure;
subplot(221)
plot(M(:,1),M(:,3),'k.','Color',[.7 .7 .7]);
hold on
plot(M(:,2),M(:,4),'k.','Color',[.7 .7 .7]);
ylim([4000 8000]);
% lsline
xlabel('Concha Size');
ylabel('Low Frequency (Hz)');
axis square
ax = axis;
x = [M(:,1);M(:,2)];
x = [x ones(size(x))];
y = [M(:,3);M(:,4)];
[b,bint,r,rint,stats] = regress(y,x);
str              = ['P = ' num2str(stats(3),2) ];
text((ax(2)-ax(1))*0.1+ax(1),(ax(4)-ax(3))*0.1+ax(3),str,'HorizontalAlignment','center');
plot([min(x(:,1));max(x(:,1))],[min(x);max(x)]*b,'k-','LineWidth',2)

subplot(222)
plot(M(:,1),M(:,5),'k.','Color',[.7 .7 .7]);
hold on
plot(M(:,2),M(:,6),'k.','Color',[.7 .7 .7]);
ylim([6000 13500]);
% lsline
xlabel('Concha Size');
ylabel('High Frequency (Hz)');
axis square
ax = axis;
x = [M(:,1);M(:,2)];
x = [x ones(size(x))];
y = [M(:,5);M(:,6)];
[b,bint,r,rint,stats] = regress(y,x);
str              = ['P = ' num2str(stats(3),2) ];
text((ax(2)-ax(1))*0.1+ax(1),(ax(4)-ax(3))*0.1+ax(3),str,'HorizontalAlignment','center');
plot([min(x(:,1));max(x(:,1))],[min(x);max(x)]*b,'k-','LineWidth',2)
set(gca,'YAxisLocation','right');

subplot(223)
x = [age;age];
y = [M(:,1);M(:,2)];
plot(x,y,'k.','Color',[.7 .7 .7]);
hold on
xlabel('Age (yrs)');
ylabel('Concha Size');
axis square
ax = axis;
x = [age;age];
x = [x ones(size(x))];
y = [M(:,1);M(:,2)];
[b,bint,r,rint,stats] = regress(y,x); %#ok<*ASGLU>
str              = ['R^2 = ' num2str(stats(1)^2,2) ', P = ' num2str(stats(3),2) ];
text((ax(2)-ax(1))*0.1+ax(1),(ax(4)-ax(3))*0.9+ax(3),str);
plot(ax([1 2]),[ax([1 2]); 1 1]'*b,'k-','LineWidth',2)

% subplot(224)
% x = [age;age;age1;age2;age3;age4;age5];
% y = [M(:,3);M(:,4);Flow1;Flow2;Flow3;Flow4;Flow5];
% plot(x,y,'k.','Color',[.7 .7 .7]);
% hold on
% xlabel('Age (yrs)');
% ylabel('Low Frequency (Hz)');
% axis square
% ax = axis;
% x = [x ones(size(x))];
% [b,bint,r,rint,stats] = regress(y,x);
% str              = ['R^2 = ' num2str(stats(1)^2,2) ', P = ' num2str(stats(3),2) ];
% text((ax(2)-ax(1))*0.1+ax(1),(ax(4)-ax(3))*0.9+ax(3),str);
% plot(ax([1 2]),[ax([1 2]); 1 1]'*b,'k-','LineWidth',2)
% set(gca,'YAxisLocation','right');

% 
% x = [age1;age2;age3;age4;age5];
% y = [Flow1;Flow2;Flow3;Flow4;Flow5];
% plot(x,y,'ko','MarkerFaceColor','k');
% 
% figure
% x = [age;age;age1;age2;age3];
% x = round(x/2)*2;
% y = [M(:,3);M(:,4);Flow1;Flow2;Flow3];
% SPKbubblegraph(x,y);