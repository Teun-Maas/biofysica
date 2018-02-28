function cfg = tdt_globals(cfg)
% TDT_GLOBALS
%
% "Globals" for TDT system:
% - Circuit filenames for RA16 and RP2s
%
% These circuit filenames are searched on the Matlab path.

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com

%% TDT flags
cfg.statTDTConnect = 1;
cfg.statTDTLoad    = 2;
cfg.statTDTRun     = 3;


%% Fixed RA16.rcx
cfg.medusaFs	= 6103.515625; % Hz
cfg.dataidx		= {'Data_1' 'Data_2' 'Data_3'}; % names of Data sources

%% RP2 and Muxes
% cfg.mux2rp2				= [1 2 2 1]; % which RP2 channel belongs to which MUX?
cfg.mux2rp2				= [1 2]; % which RP2 channel belongs to which MUX?

% cfg.recdataidx		= {'recData_1' 'recData_2'}; % names of Data sources

%% Fixed setup (SA1 & PA5 & RP2.rcx)
cfg.maxsndlevel		= 75; % ?

%% Hardware "defaults" - needs to be checked every time
cfg.SA1gain			= 0;
% 
% cfg.remmel.X.gain	= 500;
% cfg.remmel.X.offset	= 500;
% cfg.remmel.X.invert	= 0;
% cfg.remmel.Y.gain	= 500;
% cfg.remmel.Y.offset	= 500;
% cfg.remmel.Y.invert	= 0;
% cfg.remmel.Z.gain	= 500;
% cfg.remmel.Z.offset	= 500;
% cfg.remmel.Z.invert	= 0;

%% standard variables
if ~isfield(cfg,'RZ6_1circuit')
	% RP2_2circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\hrtf_2Channels_rp2.rco';
	% RP2_2circuit = 'C:\Gitlab\Marc\experiment\RPvdsEx\GvB_GWN_Single_Sound.rcx';
	cfg.RZ6_1circuit = which('sphereMinor_RZ6.rcx');
end


%% short names of TDT circuits
[~,cfg.sRZ6_1circuit]	= fileparts(cfg.RZ6_1circuit);

