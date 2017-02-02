function ddat = pa_nirs_detrend(dat)
% NIRS = PA_NIRS_DETREND(NIRS,FDOWN);
%
% Detrend DAT row-by-row, taking care of NaNs
%
% See also PA_NIRS_READ, PA_NIRS_PREPROCESS

% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

%% Initialization


%% detrend
% sel				= time<10; % remove everything before 8 s, this is based on observing the data
% oxychan(:,sel)	= NaN;

m			= nanmean(dat,2);
dat	= bsxfun(@minus,dat,m);
% mu			= nanmean(dat);
% % mu = detrend(mu);
nchan = size(dat,1);
ddat = dat; % detrended dat
for chanIdx = 1:nchan
	x = dat(chanIdx,:);
	sel			= ~isnan(x);
	if sum(sel)
	x(sel)		= detrend(x(sel));
	ddat(chanIdx,:) = x;
	end
end
