function nirs = pa_nirs_sci(nirs,varargin)
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
