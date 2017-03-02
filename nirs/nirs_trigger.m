function data = nirs_trigger(data,varargin)
% TRG = NIRS_TRIGGER(DATA,CHAN)
%
% Determine trigger index TRG for Oxymon AD channels CHAN in DATA
%
% NIRS_TRIGGER assumes that the trigger has a duration TRGDUR, and that the
% trigger has a voltage > 2.5V.
% TRG is a cell array of numel(CHAN)x1
%
% method
% chan
% fsdown
% 
%
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

meth		= keyval('method',varargin,'kennan');
fsdown		= keyval('fsdown',varargin,10);
trgdur		= keyval('trgdur',varargin,100);
chan		= keyval('chan',varargin,[97 98]);

nchan		= numel(chan);
fs			= data.fsample;

%%
if strcmpi(meth,'kennan')
trigger = cell(nchan,1); % samples
	
	for chnIdx = 1:nchan
		AD		= data.trial{1}(chan(chnIdx),:);
		trigger(chnIdx) = {gettrig(AD,trgdur,fs)};
		triggeroff(chnIdx) = {gettrig(AD,trgdur,fs)+0.1*fs};
		synch = {};
	end
elseif strcmpi(meth,'ciav')
	AD				= data.trial{1}(chan,:);
	[onset,offset,sync]	= gettrig_ciav(AD);
	trigger(1)		= {onset};
	synch(1)			= {transpose(double(sync))};
	triggeroff(1)		= {offset};
end

%% offset and "downsample"
nchan		= numel(trigger);
for ii = 1:nchan
    trigger{ii}				= trigger{ii}/fs*fsdown;
    triggeroff{ii}				= triggeroff{ii}/fs*fsdown;
end


% synch{1}				= resample(synch{1},fsdown,data.fsample);

%% Save
data.trigger	= trigger;
data.triggeroff	= triggeroff;
data.sync		= synch;
for ii = 1:nchan
    data.triggertime{ii}	= trigger{ii}/fsdown; % (s)
    data.triggerofftime{ii}	= triggeroff{ii}/fsdown; % (s)
end
data.triggertime = data.triggertime';
data.triggerofftime = data.triggerofftime';

function trigger = gettrig(AD,trgdur,fs)
% trigger for 'Kennan' experiment
% sel			= AD>2.5;
v			= [0 diff(AD)];
trigger		= find(v>0.0001);
dtrigger	= [max(trigger) diff(trigger)];
sel			= dtrigger<trgdur/1000*fs;
while sum(sel)
	dtrigger	= dtrigger(~sel);
	trigger		= trigger(~sel);
	sel			= dtrigger<trgdur/1000*fs;
end

% trigger = trigger-0.03*fs;

function [onset,offset,sync] = gettrig_ciav(sync)

% 1 = movie
% 2 = movie
% 3 = user interface
% 4 = nothing (screen refresh)
% remove glitches from screen refresh
b			= regstats(sync(1,:),sync(2,:),'linear','r'); 
sync		= b.r;
sync		= sync-sync(1); % offset

% then 'digitize' response
threshold	= 2*std(sync); % arbitrary, should be low enough to include onset, but large enough not to include glitches
sync		= sync>threshold; 

% Detect on- and offset
% use differences
d		= [0; diff(sync)];
onset	= find(d>0);
offset	= find(d<0);

% How many?
nmovie = numel(onset); % should be 120
nstim = 120;
% if nmovie~=nstim
% 	error('NMovieDoesNotMatchNStim');
% else
	disp([num2str(nstim) ' movies were presented']);
% end

