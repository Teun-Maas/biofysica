function data = pa_nirs_matread(fname)
% PA_NIRS_MATREAD reads NIRS time signals from mat-file
% format with NIRS-structure.
%
% NIRS = PA_NIRS_MATREAD(FILENAME)
%
% The output spike structure contains
%   NIRS.label		= Nchansx1 cell-array, with channel labels
%   NIRS.event		= 1xNevent struct, with stimulus events
%   NIRS.hdr		= header structure/information
%   NIRS.time		= 1xNsamples array, each element contains sample time
%   NIRS.trial      = NchansxNsamples array, each element contains a
%						measurement
%	NIRS.cfg		= configuration structure
%
% PA_NIRS_READ determines 'events' from AD channels. To do: read events
% from event channel

% (c) 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Check file
if nargin<1
	fname = pa_fcheckexist([]);
end
%% Load data from Excelsheet
nirs	= load(fname);
data	= nirs.nirs;
% [chan,TXT,RAW]	= xlsread(fname);


%% Get events
event			= read_event(data.ADvalues,data.Fs);
data.event		= event;

%% Bookkeeping
% add version information to the configuration
cfg.version.name	= mfilename('fullpath');

% add information about the Matlab version used to the configuration
cfg.version.matlab	= version;

% remember the configuration details of the input data
cfg.previous = [];
try %#ok<TRYNC>
	cfg.previous	= cfg;
end
% remember the exact configuration details in the output
data.cfg = cfg;


function event = read_event(eventchans,Fs)
% [TRIAL,EVENT] = READ_EVENT(FILENAME)
%
% READ_EVENT reads all events from a NIRS dataset and returns
% them in a well defined structure.ADindx = find(strncmpi('A/D',hdr.label,3));

% ADindx			= strncmpi('A/D',hdr.label,3);
% eventchans 			= chan(hdr.startdata+1:end,ADindx);
% eventchans = eventchans(3:end,:);
eventchans = eventchans(:,1:2); % This only works for one particular data set!
switch size(eventchans,2)
	case 1
		trial			= eventchans;
		rnd             = 0;
	case 2
		trial			= sum(eventchans,2);
		rnd             = 1;
		trial(trial>2)  = 4; % set low-above-0 values to high
	otherwise
		MsgID = 'PandA:pa_nirs_read:read_event:NumberofADChansWrong';
		ErrMsg = 'Number of AD channels is larger than 2. Please remove non-informative channels from excel-file';
		warning(MsgID,ErrMsg);
		trial			= sum(eventchans,2);
		rnd             = 1;
end
peak            = peak_detect(trial, Fs);
event			= struct([]);
for ii = 1:length(peak)
	event(ii).stim = [];
	if mod(ii,2)
		event(ii).value		= 1; % Assumption: first peak is onset. Is this correct?
	else
		event(ii).value		= 0;
	end
	event(ii).type		= 'backpanel ADC1';
	event(ii).sample	= peak(ii); % (samples), the first sample of a recording is 1, first event set at + 10 s
	event(ii).offset	= 0; % (samples)
	if ii<length(peak)
		event(ii).duration	= peak(ii+1)-peak(ii); % (samples)
	else % last offset of pulse
		event(ii).duration	= 0; % (samples)
	end
	
	if rnd && any(eventchans(peak(ii)-2:peak(ii)+2,1) > 2)
		event(ii).stim = 'A';
	end
	if rnd && any(eventchans(peak(ii)-2:peak(ii)+2,2) > 2)
		event(ii).stim = [event(ii).stim 'V'];
	end
end

function peak = peak_detect(trial, Fs)
pulse			= [0; diff(trial)];	% 'differentiate' pulse
lvl				= 4*std(pulse);		% threshold set at 4 SD
peak			= find(pulse>lvl | pulse<-lvl);	% Obtain peaks (sample numbers) both onset and offset as double check
mxvol           = floor(max(trial));

% Check to see whether peaks are not too close together
% e.g. each peak should be minimally half the average duration away
muduration	= 5;
sel			= false;
while any(~sel)
	peak1		= peak(1:end-1);
	sel			= (peak-[0;peak1])>(muduration/2);
	sel(1)      = 1; % First peak doesn't need to be > (muduration/2)
	if ~isempty(peak)
		peak		= peak(sel);
	end
end

% Check if event lasts at least 1 sec
minduration = Fs;
del = [];
for ii = 1:numel(peak)
    indx = peak(ii)-minduration:peak(ii);
	pre  = sum(trial(indx));
	pre  = Fs*round(pre/Fs);

    indx = peak(ii):peak(ii)+minduration;
sel = indx<size(trial,1);
indx = indx(sel);
    post = sum(trial(indx));
	post = Fs*round(post/Fs);
	
	tpre  = pre < minduration * mxvol;
	tpost = post < minduration * mxvol;
	if tpre && tpost || tpost && ii == 1 || tpre && ii == numel(peak)
		del = [del ii];
	end
end
peak(del) = [];