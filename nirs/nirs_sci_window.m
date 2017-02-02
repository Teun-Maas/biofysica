function data = nirs_sci_window(data,varargin)
% SCI = NIRS_SCI(DATA)
%
% Determine scalp coupling index
%
% See also Pollonini et al 2014:
% Pollonini, L., Olds, C., Abaya, H., Bortfeld, H., Beauchamp, M. S., & Oghalai, J. S. (2014).
% Auditory cortex activation to natural speech and simulated cochlear implant speech measured
% with functional near-infrared spectroscopy. Hearing Research, 309, 84?93. doi:10.1016/j.heares.2013.11.007

%% Initialization
dispFlag = keyval('disp',varargin,false);

%%
% read raw optical density
OD1		= data.trial{1};
sel		= regexp(data.label,'R*T*','once'); % search for all Receivers and Transmitters
sel		= ~cellfun(@isempty,sel);
% sum
OD1		= OD1(sel,:);
%%
% hdr = nirs.xmlInfo;
% yields 16 rows: Laser 1 to 8 for 2 detectors
Fs = data.fsample;
nChan = size(OD1,1);

%% Get cardiac component
OD			= OD1;
for ii = 1:nChan
	OD(ii,:) = resample(OD1(ii,:),10,Fs); %#ok<*AGROW,*SAGROW> % we resample the data: this is better than downsample because it deals with anti-aliasing, but there is a discussion about this, see also auditory mailing list
	OD(ii,:) = bandpass(OD(ii,:),[0.5 2.5],5); % we bandpass between 0.5 and 2.5 Hz to keep cardiac component only.
end

%% Coorelate cardiac component for two wavelengths of one optode
% The idea is that a well-placed optode should have a cardiac component
% at both wavelengths

%%
kk = 0;
SCI = NaN(round(nChan/2),1);
sb = ceil(sqrt(nChan/2));
for ii = 1:2:nChan % Optical density, so pairs of lasers
	
	kk		= kk+1; % counter
	
	% remove some outliers, which may also include filtering artefacts
	t		= prctile(OD(ii,:),[2.5 97.5]);
	sel1	= OD(ii,:)<t(2) & OD(ii,:)>t(1);
	t		= prctile(OD(ii+1,:),[2.5 97.5]);
	sel2	= OD(ii+1,:)<t(2) & OD(ii+1,:)>t(1);
	sel		= sel1 & sel2;
	
	%% zscore - basically useless?
% 	OD(ii,sel)		= zscore(OD(ii,sel)); % zscore to remove amplitude differences
% 	OD(ii+1,sel)	= zscore(OD(ii+1,sel)); % zscore to remove amplitude differences
	% cross-correlation with lag0
	% 		try %#ok<TRYNC> % yes, poor programming
	a = OD(ii,sel);
	b = OD(ii+1,sel);
	
	
avg = movcorr(a,b,100);
whos a b avg

figure(1)
clf
subplot(411)
plot(a,'.');

subplot(412)
plot(b,'.');

subplot(414)
plot(a,b,'.');
axis square
subplot(413)
plot(avg,'.')
ylim([-1.1 1.1]);
drawnow
pause
	% 			end
end
s			= [SCI SCI]'; % double for all combinations of lasers
data.sci	= s(:);

if dispFlag
	figure(3)
	hist(data.sci,0:0.05:1);
	axis square
	box off
	verline(0.75,'r-');
	xlabel('Scalp coupling index');
	ylabel('Number of channels');
end
