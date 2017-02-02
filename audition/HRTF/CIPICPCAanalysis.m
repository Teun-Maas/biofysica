8clear all
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
HL = [];
AZ = [];
GWN = randn(2^12,1);
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
	hl = squeeze(hrir_l(:,Ie,:))';               % Normalize data [-1 1]
	whos hl
	
	
	
% 	remove phase
	x = hl(:,1);
	for jj = 1:LA
		
		figure(3)
		plot(hl(:,jj)+jj);
		hold on
		[c,l] = xcorr(x,hl(:,jj),20,'coeff');
		[~,idx] = max(c);
		
% 		if jj>1
		hl(:,jj) = circshift(hl(:,jj),idx);
% 				hl(:,jj) = hl(:,jj)./max(hl(:,jj));

% 		end
% 		figure(2)
% 		clf
% 		plot(l,c)
% 		hold on
% 		title(A(jj));
% 	drawnow
% 	pause(.5)
	end
	HL = [HL hl];
	
	AZ = [AZ A];
	figure(4)
	clf
	plot(hl,'Color',[.7 .7 .7]);
	axis square;
	box off
	hold on
	plot(mean(hl,2),'k-','LineWidth',2);
	
	
	[U,S,V] = svd(hl);
	s = diag(S);
	figure(5)
	imagesc(S);
	
	
	figure(6)
	clf
	subplot(221)
	plot(A,V(:,1).*s(1));
	box off;
	axis square;
	
	subplot(222)
	plot(A,V(:,2)*s(2));
	box off;
	axis square;
	
	subplot(223)
	plot(s./sum(s),'ko');
	box off
	axis square;
% 	return
% 	figure(1)
% 	clf
% 	plot(hl)
% 	drawnow
% 	pause(.5)
end
%%


	[U,S,V] = svd(HL);
	s = diag(S);
	figure(7)
	imagesc(S);
	
	
	figure(8)
	clf
	subplot(221)
	plot(AZ,V(:,1),'k.');
	box off;
	axis square;
	
	subplot(222)
	plot(AZ,V(:,2),'k.');
	box off;
	axis square;
	
	subplot(223)
	plot(s./sum(s),'ko');
	box off
	axis square;
%%
keyboard
% only azimuth? what about elevation? is there no confound?
% show HRIRs
% what about ITD? Lag?
% m=360 from TU Berlin, which ones? description of database would be
% preferred: one subjects, one kemar-head, many subjects?
% why distance? Is not supposed to be different, except for power...
% why HSE? should you not remove this ambiguous cue (similar to the phase, or binaural)
% abbreviations, missing or unnecessary
% figure 1: is right azimuth indicated by negative values?
% decay rate of eigenvalues, how well does the 1st-component model fit the
% data?
% why only left HRIRs? Can't you pool the data (by reversing the azimuth angles)? Or can you predict right
% HRIRs?
