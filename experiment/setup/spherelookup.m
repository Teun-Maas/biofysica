function cfg = spherelookup(cfg)
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
cfg					= sphere_speakerpositions(cfg);
cfg.nstimchan		= 2^7; % number of PLC and MUX channels
cfg.nspeakers		= 112; % actual number of speakers
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
cfg.nRP2			= 2; % 2 RP2s
cfg.nMUX			= 4; % 4 multiplexers
cfg.nMUXbit			= 16; % 16 bits per MUX
cfg.RP2ind			= 1:cfg.nRP2;
cfg.MUXind			= 1:cfg.nMUX;
cfg.MUXbit			= 1:cfg.nMUXbit;
LookUp(:,1)			= cfg.stimchan; %
for ii				= 1:cfg.nRP2
	idx				= (1:cfg.nMUX*cfg.nMUXbit)+(ii-1)*cfg.nMUX*cfg.nMUXbit;
	LookUp(idx,2)	= cfg.RP2ind(ii);
	for jj			= 1:cfg.nMUX
		idx				= (1:cfg.nMUXbit)+(jj-1)*cfg.nMUXbit+(ii-1)*cfg.nMUX*cfg.nMUXbit;
		LookUp(idx,3)	= cfg.MUXind(jj);
	end
end
LookUp(:,4)			= repmat(cfg.MUXbit,1,cfg.nRP2*cfg.nMUX);

%% Combine
% [PLC-Chan# RP2# MUX# BIT# AZ EL]
cfg.lookup		= [cfg.lookup(:,3) LookUp(:,2:4) cfg.lookup(:,1)  cfg.lookup(:,2)];
cfg.lookuplabel = {'Channel' 'RP2' 'MUX' 'Bit' 'Azimuth' 'Elevation'};

%% Channel 79 is missing
% quick fix, replace with Channel 31
% Date: 20-8-2015
cfg.lookup(32,5:6) = cfg.lookup(80,5:6);
cfg.lookup(80,5:6) = NaN;

%% Extra infrared
cfg.lookup(128,5:6) = [700,700];
% keyboard
%% Tannoy speakers
% Januray 8
TannoyLocations				= [-70:10:-10 10:10:70];
TannoyIDs					= [26:30 58:63  124:126];
cfg.TannoyIDs				= TannoyIDs;

cfg.lookup(TannoyIDs+1,5)	= TannoyLocations; % azimuth
cfg.lookup(TannoyIDs+1,6)	= +7.5; % elevation


%% Interpolant
x		= cfg.lookup(:,5);
y		= cfg.lookup(:,6);
z		= cfg.lookup(:,1);
sel		= ~isnan(x);
x		= x(sel);
y		= y(sel);
z		= z(sel);
cfg.interpolant		= scatteredInterpolant(x,y,z,'nearest');