function cfg = readcfg(cfg)
% CFG = READCFG(CFG)
%
% Load ConFiGuration parameters from exp-file FNAME
%
% See also READEXP

% 2015 Marc van Wanrooij

%% Initialization
if nargin<1
	cfg.cfgfname = which('HumanV1.cfg');
% 	cfg		= getexp;
% 	fname	= cfg.expfname;
end
if isstruct(cfg)
	fname = [cfg.expdir filesep cfg.cfgfname];
end
%% Open
fid		= fopen(fname,'r');

%% Skip first three comment lines
for trlIdx = 1:5
	fgetl(fid);
end
fscanf(fid,'%s',1); % Data folder
cfg.humanv1.datmap = fscanf(fid,'%s',1);
cfg.humanv1.datmap = cfg.humanv1.datmap(2:end-1);

fscanf(fid,'%s',1); % EXP file folder
cfg.humanv1.expmap = fscanf(fid,'%s',1);
cfg.humanv1.expmap = cfg.humanv1.expmap(2:end-1);

fscanf(fid,'%s',1); % sound file folder
cfg.humanv1.sndmap = fscanf(fid,'%s',1);
cfg.humanv1.sndmap = cfg.humanv1.sndmap(2:end-1);

for trlIdx = 1:3
	fgetl(fid);
end

%% TDT
fscanf(fid,'%s',1); % RP2 circuit
cfg.humanv1.RP2_1circuit = fscanf(fid,'%s',1);
cfg.humanv1.RP2_1circuit = cfg.humanv1.RP2_1circuit(2:end-1);

fscanf(fid,'%s',1); % RP2 circuit
cfg.humanv1.RP2_2circuit = fscanf(fid,'%s',1);
cfg.humanv1.RP2_2circuit = cfg.humanv1.RP2_2circuit(2:end-1);

fscanf(fid,'%s',1); % RA16 circuit
cfg.humanv1.RA16circuit = fscanf(fid,'%s',1);
cfg.humanv1.RA16circuit = cfg.humanv1.RA16circuit(2:end-1);

for trlIdx = 1:3
	fgetl(fid);
end

%% Channels
for ii = 1:8
cfg.humanv1.ADC(ii).channel		= fscanf(fid,'%s',1);
cfg.humanv1.ADC(ii).lp			= fscanf(fid,'%d',1);
cfg.humanv1.ADC(ii).rate		= fscanf(fid,'%d',1);
cfg.humanv1.ADC(ii).samples		= fscanf(fid,'%d',1);
cfg.humanv1.ADC(ii).name		= fgetl(fid);
end

%% Close
fclose(fid);
