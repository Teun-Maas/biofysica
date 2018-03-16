function cfg = spherelookupMinor(cfg)
% CFG = SPHERELOOKUP
%
% Get lookup table for LEDs and speakers in the Sphere setup
%
% CFG contains field
%	- lookup: % [PLC-Chan# RP2# MUX# BIT# AZ EL]

if nargin<1
	cfg = [];
end
%% LED connections
% [Azimuth (deg) Elevation (deg) LED #]
% azimuths and elevations were measured by Sebastian Ausili in June 2015
cfg					= sphere_speakerpositionsMinor(cfg);
cfg.nstimchan		= 2^5; % number of PLC and MUX channels
cfg.nspeakers		= 29; % actual number of speakers
cfg.stimchan		= (1:cfg.nstimchan)-1;

%% Add missing channels
n					= size(cfg.lookup,1);
cfg.lookup			= [cfg.lookup; NaN(cfg.nstimchan-n,3)]; % add missing channel data as NaNs
sel					= ismember(cfg.stimchan,cfg.lookup(:,3)); % lookup existing channels
cfg.missingchan		= cfg.stimchan(~sel); % get missing channels
sel					= isnan(cfg.lookup(:,3)); % search for missing channels in lookup-table
cfg.lookup(sel,3)	= cfg.missingchan; % put missing channel numbers in lookup-table
cfg.lookup			= sortrows(cfg.lookup,[3 1 2]); % sort lookup-table by channel number

%% Speakers
% [Channel# RP2# MUX# BIT#]
%
% RP2ind   = [1 2];
% MUXind   = 1:4;
% MUXbit   = 1:16;
% ChanNo  = 0:127;
% LookUp = NaN(128.4);
% LookUp(:.1) = ChanNo;
% LookUp(1:64.2) = RP2ind(1);
% LookUp(65:128.2) = RP2ind(2);
% LookUp(1:16.3) = MUXind(1);
% LookUp(17:32.3) = MUXind(2);
% LookUp(33:48.3) = MUXind(3);
% LookUp(49:64.3) = MUXind(4);
% LookUp(65:80.3) = MUXind(1);
% LookUp(81:96.3) = MUXind(2);
% LookUp(97:112.3) = MUXind(3);
% LookUp(113:128.3) = MUXind(4);
% LookUp(:.4) = [MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit MUXbit]';

% This is still written to RP2 setup
% cfg.nRP2			= 2; % 2 RP2s
% cfg.nMUX			= 4; % 4 multiplexers
% cfg.nMUXbit			= 16; % 16 bits per MUX
% cfg.RP2ind			= 1:cfg.nRP2;
% cfg.MUXind			= 1:cfg.nMUX;
% cfg.MUXbit			= 1:cfg.nMUXbit;

cfg.nMUX			= 2; % 2 multiplexers
cfg.nMUXbit			= 16; % 16 bits per MUX
cfg.MUXind			= 1:cfg.nMUX;
cfg.MUXbit			= 1:cfg.nMUXbit;

LookUp(:,1)			= cfg.stimchan; % Stimulus channel 0-31
idx					= LookUp(:,1)+1;
LookUp(:,2)			= ceil(idx/(2^4));% MUX number
idx					= mod(idx,2^4);
idx(idx==0)			= 2^4;
LookUp(:,3)			= idx; % TDT bit

%% Combine
% [PLC-Chan# RP2# MUX# BIT# AZ EL]
cfg.lookup		= [cfg.lookup(:,3) LookUp(:,2:3) cfg.lookup(:,1)  cfg.lookup(:,2)];
cfg.lookuplabel = {'Channel' 'MUX' 'Bit' 'Azimuth' 'Elevation'};
% cfg.lookuplabel = {'Channel' 'RP2' 'MUX' 'Bit' 'Azimuth' 'Elevation'};


%% Interpolant
x		= cfg.lookup(:,4);
y		= cfg.lookup(:,5);
z		= cfg.lookup(:,1);
sel		= ~isnan(x);
x		= x(sel);
y		= y(sel);
z		= z(sel);
cfg.interpolant		= scatteredInterpolant(x,y,z,'nearest');