function cfg = tdt_globalsMinor(cfg)
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


%% Fixed RZ6.rcx
% cfg.RZ6Fs	= 48828.125; % Hz = default, will be replaced in 
cfg.dataidx		= {'Data_1' 'Data_2' 'Data_3'}; % names of Data sources

%% RP2 and Muxes
% cfg.mux2rp2				= [1 2 2 1]; % which RP2 channel belongs to which MUX?
cfg.mux2rp2				= [1 2]; % which RP2 channel belongs to which MUX?

% cfg.recdataidx		= {'recData_1' 'recData_2'}; % names of Data sources

%% Fixed setup (SA1 & PA5 & RP2.rcx)
% cfg.maxsndlevel		= 75; % ?

%% standard variables
if ~isfield(cfg,'RZ6_1circuit')
	% RP2_2circuit = 'C:\matlab\Marc\Experiment\RPvdsEx\hrtf_2Channels_rp2.rco';
	% RP2_2circuit = 'C:\Gitlab\Marc\experiment\RPvdsEx\GvB_GWN_Single_Sound.rcx';
	cfg.RZ6_1circuit = which('sphereMinor_RZ6_mp.rcx');
end


%% short names of TDT circuits
[~,cfg.sRZ6_1circuit]	= fileparts(cfg.RZ6_1circuit);

