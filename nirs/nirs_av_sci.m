function nirs_av_sci
close all


% scalpFlag = true;
scalpFlag = false;

%% Scalp coupling
if scalpFlag
	k= 0;
	% Load mat files
	cd('/Users/marcw/DATA/NIRS/OXY3_v14112014'); %#ok<*UNRCH> % contains all relevant data files
	d		= dir('0*.oxy3');
	[m,~]	= size(d);
	missingfile = cell(1);
	for ii	= 1:m
		clf
		fname	= d(ii).name;
		try  % bad programming, figure out problems later
			nirs = oxysoft2matlab(fname,'rawOD',[]);
			nirs = pa_nirs_sci(nirs,'disp',false);
			drawnow
			
			fname = pa_fcheckext(fname,'mat');
			save(fname,'nirs')
		catch
			k=k+1;
			missingfile{k} = fname;
		end
	end
end


%% Overview xlsx
% %% Scalp coupling
% % Load mat files
% cd('/Users/marcw/DATA/NIRS/OXY3_v14112014'); %#ok<*UNRCH> % contains all relevant data files
% d		= dir('0*.xlsx');
% [m,~]	= size(d);
% % SCI = struct([]);
% for ii	= 1:m
% 	clf
% 	fname	= d(ii).name;
% 	disp(fname);
% 	nirs = pa_nirs_readxls(fname);
% 	nirs
% 	% 		try %#ok<TRYNC> % bad programming, figure out problems later
% 	% 			nirs = oxysoft2matlab(fname,'rawOD',[]);
% 	% 			nirs = pa_nirs_sci(nirs,'disp','false');
% 	% 			drawnow
% 	%
% 	% 			fname = pa_fcheckext(fname,'mat');
% 	% 			save(fname,'nirs')
% 	% 		end
% end
% return

%% Preprocessing
cd('/Users/marcw/DATA/NIRS/OXY3_v14112014'); %#ok<*UNRCH> % contains all relevant data files
d		= dir('0*.mat');

[m,~]	= size(d);
SCI = [];
F = [];
nsci = NaN(m,1);
for ii	= 1:m
	% 		for ii = 39
	fname		= d(ii).name;
	disp(fname)
	nirs		= pa_nirs_matread(fname); % read in mat file, and detect events/stimuli/btn from AD channels
	sci			= nirs.sci;
	nsci(ii)	= numel(sci);
	SCI			= [SCI;sci];
	F			= [F;repmat(str2num(fname(1:3)),nsci(ii),1)]; %#ok<ST2NM>
end


%% Overview oxy3
median(SCI)

[m numel(SCI) sum(SCI<=0.9)]



%% Check
cd('/Users/marcw/DATA/NIRS/OXY3_v14112014');

d = dir('0*.oxy3');
numel(d)
l = char({d.name}');
prel = l(:,1:3);
ul = unique(prel,'rows');
n = size(ul,1);
np = size(prel,1);

K = NaN(n,1);
for ii = 1:n
	k = 0;
	for jj = 1:np
		t = strfind(prel(jj,:),ul(ii,:));
		if ~isempty(t)
		k = k+1;
		end
	end
	K(ii) = k;
end

nside = (str2num(ul)<38 & str2num(ul)>15)+1; %#ok<ST2NM>
nchan = sum(K.*nside*2);% Ref vs Sig & number sides
disp(nchan)

% %% Overview xlsx
% %% Scalp coupling
% 	% Load mat files
% 	cd('/Users/marcw/DATA/NIRS/OXY3_v14112014'); %#ok<*UNRCH> % contains all relevant data files
% 	d		= dir('0*.xlsx');
% 	[m,~]	= size(d);
% 	SCI = struct([]);
% 	for ii	= 1:m
% 		clf
% 		fname	= d(ii).name
% 		try %#ok<TRYNC> % bad programming, figure out problems later
% 			nirs = oxysoft2matlab(fname,'rawOD',[]);
% 			nirs = pa_nirs_sci(nirs,'disp','false');
% 			drawnow
% 			
% 			fname = pa_fcheckext(fname,'mat');
% 			save(fname,'nirs')
% 		end
% 	end
% 
% 
%%

a = str2num(prel);
% find(a==31); %missing
a(69) = [];

[a nsci]

% First 16 subjects were measured on one side, also 5 CI (measured by
% Hai-Yin and Louise_
% the other 21 subjects were measured on both sides (by Luuk)
% one oxy3 file (31-002) cannot be read, and has no mat file
% 3, 11, 12, 16, 26, 33-3 only have xls files

%%
nchan  = 21*2*2+21*2*1
%%
keyboard


function nirs = pa_nirs_sci2(nirs,varargin)
% SCI = PA_NIRS_SCALPCOUPLINGINDEX(NIRS)
%
% Determine scalp coupling index
%
% See also Pollonini et al 2014:
% Pollonini, L., Olds, C., Abaya, H., Bortfeld, H., Beauchamp, M. S., & Oghalai, J. S. (2014). 
% Auditory cortex activation to natural speech and simulated cochlear implant speech measured 
% with functional near-infrared spectroscopy. Hearing Research, 309, 84?93. doi:10.1016/j.heares.2013.11.007

%% Initialization
dispFlag = pa_keyval('disp',varargin);
if isempty(dispFlag)
	dispFlag = false;
end

% read raw optical density
% [OD1,hdr] = oxy3read_function; % read oxy3: bloody annoying as we cannot supply filename input...
OD1 = nirs.OD;
% hdr = nirs.xmlInfo;
% yields 16 rows: Laser 1 to 8 for 2 detectors
Fs = nirs.Fs;
nChan = size(nirs.Rx_TxId,2);

%% Get cardiac component
for ii = 1:nChan
	OD(:,ii) = resample(OD1(:,ii),10,Fs); %#ok<*AGROW,*SAGROW> % we resample the data: this is better than downsample because it deals with anti-aliasing, but there is a discussion about this, see also auditory mailing list
	OD(:,ii) = pa_bandpass(OD(:,ii),[0.5 2.5],5); % we bandpass between 0.5 and 2.5 Hz to keep cardiac component only.
end

%% Coorelate cardiac component for two wavelengths of one optode
% The idea is that a well-placed optode should have a cardiac component
% at both wavelengths
kk = 0;
SCI = NaN(round(nChan/2),1);
for ii = 1:2:nChan % Optical density, so pairs of lasers
	
	kk		= kk+1; % counter
	
	% remove some outliers, which also include filtering artefacts
	t		= prctile(OD(:,ii),[2.5 97.5]);
	sel1	= OD(:,ii)<t(2) & OD(:,ii)>t(1);
	t		= prctile(OD(:,ii+1),[2.5 97.5]);
	sel2	= OD(:,ii+1)<t(2) & OD(:,ii+1)>t(1);
	sel		= sel1 & sel2;
	
	%
	OD(sel,ii)		= zscore(OD(sel,ii)); % zscore to remove amplitude differences
	OD(sel,ii+1)	= zscore(OD(sel,ii+1)); % zscore to remove amplitude differences
	
	if dispFlag
		figure(1)
		subplot(2,4,kk)
		plot(OD1(:,ii))
		box off;
		axis square;
		ylim([0 6.5]);
		xlim([0 length(OD1(:,ii))]);
		
		% X-Y plot
		figure(2)
		subplot(2,4,kk)
		hold on
		plot(OD(sel,ii),OD(sel,ii+1),'k.');
		box off
		axis square
		ylim([-3 3]);
		xlim([-3 3]);
		pa_unityline;
		xlabel('Amplitude Laser 1 (au)');
		ylabel('Amplitude Laser 2 (au)');
	end
	% cross-correlation with lag0
	% 		try %#ok<TRYNC> % yes, poor programming
	r		= corrcoef(OD(sel,ii),OD(sel,ii+1));
	r		= r(2)^2;
	% Keep data
	SCI(kk) = r;
	% graphics
	% for setup left-right Rx1-Tx1-2, Rx2-Tx3-4
	if dispFlag
		
		% time -X/Y
		figure(3)
		subplot(2,4,kk)
		plot(OD(sel,ii+1),'k.-');
		hold on
		plot(OD(sel,ii),'r.-');
		xlim([200 300]);
		ylim([-2 2]);
		title(['SCI = ' num2str(r,2)])
		xlabel('Time (samples)');
		ylabel('Amplitude (au)');
		
		axis square;							% publication-quality
		box off;								% publication-quality
		set(gca,'TickDir','out','TickLength',[0.005 0.025],...
			'XTick',200:50:300,'YTick',-2:2,...
			'FontSize',15,'FontAngle','Italic','FontName','Helvetica'); % publication-quality
	end
	% 		end
end
nirs.sci = SCI;