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
	
	c = ['/Users/marcw/Gitlab/thirdparty/CIPIC/standard_hrir_database/subject_' num2str(sub,'%03d')];
	cd(c);
	load hrir_final
	
	
	%%
	[m,n,k] = size(hrir_l)
		hl		= squeeze(hrir_l(:,9,:))';               % Normalize data [-1 1]
			[AL]	= freq_resp(hl,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
nf = size(AL,1);
AL = NaN(nf,m,n);
	AR = NaN(nf,m,n);
	
	for Ie = 1:LE
		hl		= squeeze(hrir_l(:,Ie,:))';               % Normalize data [-1 1]
		hr		= squeeze(hrir_r(:,Ie,:))';
		[AL(:,:,Ie)]	= freq_resp(hl,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
		[AR(:,:,Ie),F]	= freq_resp(hr,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	end
% 	AL = AR;
	%%
	figure(1)
	clf
	cnt = 0;
	
	for ff = 4000:2000:10000
		cnt = cnt+1;
		[mn,idx] = min(abs(F-ff));
		closestValue = F(idx);
		
		mtx = squeeze(AL(idx,:,:));
		
		subplot(2,2,cnt)
		
		contourf(A,E,mtx');
		nicegraph;
		axis normal
		title(ff)
	end
	datadir
	savegraph(['dtfcontour' num2str(ii)],'eps');	
	
	%%
	AL = reshape(AL,nf,m*n);
	muAL = mean(AL,2);
	DTF_L = AL-muAL;
	AL = reshape(AL,nf,m,n);
	
	
	
	figure(2)
	clf
	subplot(211)
	plot(F,muAL)
	
	subplot(212)
	plot(F,DTF_L);

		datadir
	savegraph(['hr_d_tf' num2str(ii)],'eps');
	
	%%
	AL		= reshape(AL,nf,m*n);
	muAL	= mean(AL,2);
	DTF_L	= AL-muAL;
	AL		= reshape(AL,nf,m,n);
	DTF_L		= reshape(DTF_L,nf,m,n);
	
	sel=E==0;
	DTF_L_A = squeeze(DTF_L(:,:,sel));
	[mn,idx] = min(DTF_L_A');
	figure(3)
	clf
	subplot(121)
	semilogx(F,A(idx),'ko-','MarkerFaceColor','w');
	nicegraph;
	ylim([0 90]);
	set(gca,'XTick',oct2bw(1000,-1:4),'XTickLabel',oct2bw(1000,-1:4)/1000);
	xlim([250 20000]);
	xlabel('Frequency (Hz)');
	ylabel('acoustical axis (deg)');
	title('azimuth');
	text(400,80,'A');
	
	sel=A==0;
	DTF_L_E = squeeze(DTF_L(:,sel,:));
	[mn,idx] = min(DTF_L_E');
	subplot(122)
	semilogx(F,E(idx),'ko-','MarkerFaceColor','w');
	nicegraph;
	set(gca,'XTick',oct2bw(1000,-1:4),'XTickLabel',oct2bw(1000,-1:4)/1000);
	xlim([250 20000]);
	ylim([-90 90]);
	xlabel('Frequency (Hz)');
	ylabel('acoustical axis (deg)');
		title('elevation');
	text(400,80,'B');

	datadir
	savegraph(['acoaxi' num2str(ii)],'eps');
	%%
	% 	[AL]	= freq_resp(hrir_l,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	% 	[AR,F]	= freq_resp(hrir_r,Fmin*1000,Fmax*1000,Q,log_scale, NFFT, fs);
	

end
