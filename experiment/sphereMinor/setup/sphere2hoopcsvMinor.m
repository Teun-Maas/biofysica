function sphere2hoopcsv(fname)
% SPHERE2HOOPCSV(FNAME)
%
% Convert the SPHERE-file structure to CSV-file of the
% HOOP
%
% See also SPHEREMAT2HOOPDAT

if nargin<1
% 	fname = [];
	fname = 'EX-0000-15-06-23-0000.mat';
end

%% Initialization
ext             = '.sphere';
matfile         = fcheckexist(fname,ext);
matfile         = fcheckext(matfile,ext);


ext             = '.csv';
csvfile         = fcheckext(matfile,ext);
dlm				= ';';

%% Load data
load(matfile,'-mat')

%% Header

ntrials_norep	= cfg.ntrials/cfg.Repeats;

% first line, general configuration
first			= [0 ntrials_norep cfg.Repeats cfg.ntrials cfg.ITI' cfg.Random cfg.nchan];
dlmwrite(csvfile,first,dlm)

% next lines, channel configuration
M				= NaN(cfg.nchan,6);
for ii			= 1:cfg.nchan
	M(ii,:)			= [0 ii ii round(cfg.medusaFs) round(cfg.medusaFs) cfg.nsamples];
end
dlmwrite(csvfile,M,'-append','delimiter',dlm)

%% Write trial information
fid = fopen(csvfile,'a+');
for trlIdx = 1:cfg.ntrials
	nstim = trial(trlIdx).nstim;
	for stmIdx = 1:nstim
		fprintf(fid,'%d%c', trlIdx,dlm); % 1) trial number
		fprintf(fid,'%d%c',stmIdx,dlm); % 2) stimulus number in trial
		fprintf(fid,'%d%c',round(mean(cfg.ITI)),dlm); % 3) Desired Inter Trial Interval
		fprintf(fid,'%d%c',round(mean(cfg.ITI)),dlm); % 4) Actual inter trial interval
		switch trial(trlIdx).stim(stmIdx).modality
			case 'LED'
				mod = 'led';
				attr = 0; % default
			case 'sound';
				mod = 'snd1';
				attr = trial(trlIdx).stim(stmIdx).matfile;
				attr = str2double(attr(4:end-4));
			case 'trigger'
				mod = 'trg0';
				attr = 0; % default
			case 'data acquisition'
				mod = 'acq';
				attr = 0; % default
		end
		fprintf(fid,'%s%c',mod,dlm);% 5) Modality of stimulus
		fprintf(fid,'%f%c',trial(trlIdx).stim(stmIdx).azimuth,dlm);%  6) azimuth location
		fprintf(fid,'%f%c',trial(trlIdx).stim(stmIdx).elevation,dlm); % 7) elevation location
		fprintf(fid,'%d%c',trial(trlIdx).stim(stmIdx).ondelay,dlm); % 8) onset
		fprintf(fid,'%d%c',trial(trlIdx).stim(stmIdx).offdelay,dlm); % 9) offset
		fprintf(fid,'%d%c',trial(trlIdx).stim(stmIdx).intensity,dlm);% 10) intensity
		fprintf(fid,'%d%c',attr,dlm); % 11)
		fprintf(fid,'%d%c',0,dlm); % 12)
		fprintf(fid,'%d \n',1); % 13)
		
	end
end

fclose(fid);
