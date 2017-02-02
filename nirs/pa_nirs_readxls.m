function data = pa_nirs_read(fname)
% PA_NIRS_READ reads NIRS time signals from OxyMon excel data
% format.
%
% NIRS = PA_NIRS_READ(FILENAME)
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
[chan,TXT,RAW]	= xlsread(fname);

%% Load header
hdr				= read_hdr_nirs(TXT,RAW);

%% Get events
data.fsample	= hdr.Fs; % sampling frequency in Hz, single number
event			= read_event(chan,hdr,data.fsample);
data.label		= hdr.label;
data.event		= event;
data.hdr		= hdr;

%% Remove the header-information
chan			= chan((hdr.startdata+1):end,:);

%% Raw data
data.trial		= chan(:,1:hdr.nChans)'; % cell-array containing a data matrix for each trial (1 X Ntrial), each data matrix is Nchan X Nsamples
data.time		= chan(:,1)'/data.fsample; % cell-array containing a time axis for each trial (1 X Ntrial), each time axis is a 1 X Nsamples vector

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

function hdr = read_hdr_nirs(TXT,RAW)
% READ_HDR_NIRS reads header/stimulus information from NIRS time signals of
% OxyMon excel data format.
%
% HDR = READ_HDR_NIRS(FILENAME)
%
% Returns a header structure with the following elements
%   hdr.Fs                  sampling frequency
%   hdr.nChans              number of channels
%   hdr.nSamples            number of samples per trial
%   hdr.label               cell-array with labels of each channel
%
% And some additional info, which is probably never used in data analysis.
%
%   hdr.nSamplesPre         number of pre-trigger samples in each trial
%   hdr.nTrials             number of trials
%
% See also READ_HEADER, READ_NIRS

% 2011 Marc van Wanrooij

%% Insert header information in hdr
hdr.oxyexport		= TXT{1,2};
hdr.xcldate			= TXT{2,2};
hdr.filedate		= TXT{3,2};
hdr.Fs				= RAW{28,2}; %after export to xcl
hdr.samplerate		= RAW{5,2}; % (Hz)
hdr.duration		= RAW{6,2}; % (s)
hdr.nSamples		= RAW{7,2};	% number of samples per trial
hdr.optodetemplate	= RAW{10,2};
hdr.optodedistance  = RAW{11,2}; % (mm)


switch hdr.optodetemplate
	case '8 channel (split)'
		hdr.pos_indx	= [4 1 8 6 3 2 5 7];
		if hdr.optodedistance ~=40
			str = char(strcat({('--------------------------------------------------------');...
				['    According to header, Optode Distance is ' num2str(hdr.optodedistance) ' mm.'] ;...
				('    Set to 40 mm, as this makes more sense.');...
				('--------------------------------------------------------')}));
			disp(str)
			hdr.optodedistance = 40;
		end
	case '6 channel (plus 2 chan split)'
		hdr.pos_indx	= [1 5 2 6 3 7 4 8];
	case '2 chan-Sevy'
		hdr.pos_indx	= [1 2];
	case '2 x 1 channel'
		hdr.pos_indx	= [1 2];
	case '3 x 1 channel III'
		hdr.pos_indx	= [1 2 3];
	case '4 x 1 channel (split) I'
		hdr.pos_indx	= [1 2 3 4];
	otherwise
		disp('No template found');
end

% if hdr.optodedistance ~=40
% 	str = char(strcat({('--------------------------------------------------------');...
% 		['    According to header, Optode Distance is ' num2str(hdr.optodedistance) ' mm.'] ;...
% 		('    Set to 40 mm, as this makes more sense.');...
% 		('--------------------------------------------------------')}));
% 	disp(str)
% 	hdr.optodedistance = 40;
% end
hdr.DPF				= RAW{12,2}; % Differential path length
hdr.deviceid		= RAW{14,2}; % ???
hdr.nreceivers		= RAW{15,2};
hdr.nlasers			= RAW{16,2};
hdr.ADC				= RAW{17,2};
hdr.laserwavelength = [RAW{19:18+hdr.nlasers,2}];

%% Now it gets tricky:
% Have we set an export sample time and rate?
if strcmpi(TXT{29,1},'Export timespan') % Yes we do!
	hdr.exporttimespan		= [RAW{29,2:3}];
	hdr.exportsamplespan	= [RAW{30,2:3}];
end
% Start of legend varies depending on whether we set an export time
legindx = find(strcmpi('Legend',TXT(:,1)));
% We have to check what dataformat we have exported to
switch TXT{legindx+1,2}
	case 'Trace (Measurement)'
		hdr.dataformat = 'oxy'; % oxy and deoxy concentrations
	case 'Rx'
		hdr.dataformat = 'OD';	% optical densities or raw values
end
% Determine the number of data channels (columns)
leg_end_indx	= find(isnan([RAW{legindx+2:legindx+30,1}]),1);
hdr.nChans		= leg_end_indx-1;
%
switch hdr.dataformat
	case'oxy'
		for ii = 1:hdr.nChans
			hdr.label{ii,1} = RAW{legindx+1+ii,2}; % cell-array containing strings, Nchan X 1
		end
	case 'OD'
		for ii = 1:hdr.nChans
			hdr.label{ii,1} = RAW{legindx+1+ii+ii,2}; % cell-array containing strings, Nchan X 1
		end
end
% Find oxyhemoglobin
indx = NaN(hdr.nChans,1);
for ii = 1:hdr.nChans
	tmp = strfind(hdr.label{ii},'O2Hb');
	if ~isempty(tmp)
		indx(ii) = tmp;
	end
end
hdr.oxyindx		= find(~isnan(indx)); % first label is

% Find deoxyhemoglobin
indx = NaN(hdr.nChans,1);
for ii = 1:hdr.nChans
	tmp = strfind(hdr.label{ii},'HHb');
	if ~isempty(tmp)
		indx(ii) = tmp;
	end
end
hdr.deoxyindx		= find(~isnan(indx));

% Start data column
hdr.startdata = leg_end_indx+legindx+3-4;

%% Default values
hdr.nSamplesPre		= 0;
hdr.nTrials			= 1;

function event = read_event(chan,hdr, Fs)
% [TRIAL,EVENT] = READ_EVENT(FILENAME)
%
% READ_EVENT reads all events from a NIRS dataset and returns
% them in a well defined structure.ADindx = find(strncmpi('A/D',hdr.label,3));
ADindx			= strncmpi('A/D',hdr.label,3);
eventchans 			= chan(hdr.startdata+1:end,ADindx);
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
subplot(211)
plot(trial)

subplot(212)
s = repmat((1:3),size(eventchans,1),1);
whos s eventchans
plot(eventchans+s)
% title(size(eventchans,2))
peak            = peak_detect(trial, Fs);
peak
npeak = numel(peak)
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