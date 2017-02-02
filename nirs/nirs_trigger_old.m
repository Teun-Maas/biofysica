function data = nirs_trigger(data,chan,trgdur,fsdown,offset)
% TRG = NIRS_TRIGGER(DATA,CHAN,TRGDUR)
%
% Determine trigger index TRG for Oxymon AD channels CHAN in DATA
%
% NIRS_TRIGGER assumes that the trigger has a uration TRGDUR, and that the
% trigger has a voltage > 2.5V.
% TRG is a cell array of numel(CHAN)x1
%
% See also NIRS_RCS, NIRS_RMVSAMPLES

%% Initialization
if nargin<1
	w			= what('LR-04-2015-06-17_active');
	DataFolder	= w.path;
	cd(DataFolder)
	fname		= 'data_raw.mat';
	load(fname);
	data = data_raw;
end
if nargin<2
	chan = [97 98];
end
if nargin<3
	trgdur		= 100; % ms
end
if nargin<4
	fsdown = 10;
end
if nargin<5
	offset = 4000;
end

nchan		= numel(chan);
fs			= data.fsample;

%%
trigger = cell(nchan,1); % samples
for chnIdx = 1:nchan
	AD		= data.trial{1}(chan(chnIdx),:);
	trigger(chnIdx) = {gettrig(AD,trgdur,fs)};
end

%% offset and "downsample"
trigger				= cellfun(@(x) round(x/data.fsample*fsdown),trigger,'UniformOutput',false); % downsample
trigger				= cellfun(@(x) x-offset, trigger,'UniformOutput',false); % remove samples, after downsampling
triggertime			= cellfun(@(x) x/fsdown,trigger,'UniformOutput',false); % trigger (s)

%% Save
data.trigger		= trigger;
data.triggertime	= triggertime; % (s)

function trigger = gettrig(AD,trgdur,fs)
sel			= AD>2.5;
v			= [0 diff(sel)];
trigger		= find(v>0);
dtrigger	= [max(trigger) diff(trigger)];
sel			= dtrigger<trgdur/1000*fs;
while sum(sel)
	dtrigger	= dtrigger(~sel);
	trigger		= trigger(~sel);
	sel			= dtrigger<trgdur/1000*fs;
end

