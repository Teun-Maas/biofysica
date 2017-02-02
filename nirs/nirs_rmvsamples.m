function [data,nrm] = nirs_rmvsamples(data,nrm)
% DATA = NIRS_RMVSAMPLES(DATA,N)
%
% Remove first N samples from DATA.TRIAL
%
% See also NIRS_OPTODELAYOUT, NIRS_OPTODEDISTANCE, NIRS_RCS, NIRS_TRIGGER

%% Initialization

if nargin<1
	close all
	w			= what('LR-04-2015-06-17_active');
	DataFolder	= w.path;
	cd(DataFolder)
	fname		= 'data_trans.mat';
	load(fname);
	data		= data_trans;
end
if nargin<2
	nrm				= []; % number of samples to remove from start, hardcoded
end

%%
if isempty(nrm)
	h = figure(666);
	plot(data.trial{1}(1,:));
	hold on
	plot(data.trial{1}(97,:));

	[nrm,~] = ginput(1);
	nrm = round(nrm);
	close(h);
end
%% Remove samples

dat				= data.trial{:};
dat				= dat(:,nrm:end); % Hardcoded removal of first recording samples

time			= data.time{:};
time			= time(nrm:end)-nrm/data.fsample;


%% overwrite
data.time{1}	= time;
data.trial{1}	= dat;
data.nrm		= nrm;
