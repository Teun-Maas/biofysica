function slims_saccade_preprocess(fname,varargin)
% SLIMS_SACCADE_PREPROCESS(FNAME)
%
% Resamples and reshapes PupilLabs structure data to a matrix.
%
% See also: SLIMS_SACCADE_PARADIGM, RESHAPE, RESAMPLE, SACDET

% TODO:
% - convert 'parameters' to 'stimuli','configuration','data'
% - remove poor data (or incorporate into sacdet)
% - targets in correct coordinate system

%% Initialization
if nargin<1
	cd('/Users/marcw/Dropbox/Manuscript/stage/2122/Jesse Janssen/data/001-22-08-03');
	% 	fname = 'JJ-001-004-06-17-22.mat';
	% 	fname = 'JJ-001-005-06-17-22.mat';
	% 	fname = 'JJ-001-006-06-17-22.mat';
	fname = 'JJ-001-22-08-03-000.mat';
end

dsp = keyval('display',varargin,true);

%% Load
% load(fname,'stimuli','configuration','data');
load(fname,'data','configuration','stimuli');

pldata			= data.pupildata; % what parameters are in here?
evdata			= data.eventdata;

tar				= stimuli.target;

%% get time stamps and align
tpl				= lsl_correct_lsl_timestamps(pldata);
tev				= lsl_correct_lsl_timestamps(evdata);

% tpl = pldata.Timestamps;
% tev = evdata.Timestamps;

ntrials			= length(tev)/2;
fs				= 200;

%% at the moment the stamps for the fixation and the target are all in one vector
% remove the time offset (1st sample)
% t0				= tpl(1);
t0				= tev(1);


tpl				= tpl - t0;
events			= tev - t0;
tareve			= events(2:2:end); % event start for the 2nd dot = target in saccade paradigm
fixeve			= events(1:2:end); % event start for the 1st dot = central fixation in saccade paradigm

% What is each column???
HVF				= pldata.Data(13:15,:)'; % horizontal, vertical and frontal?
pl_confidence	= pldata.Data(1,:); % confidence levels

if dsp
	figure(42)
	clf
	plot(tpl,HVF,'linewidth',2);
	hold on

	plot(tpl,pl_confidence,'.','linewidth',2);
	ylim([-1.1 1.1]);
	verline(tareve,'k--');
	verline(fixeve,'g--');
	title('subject 002, session 005, 17-06-2022,');
	legend('horizontal','vertical','frontal','confidence','fixation_{onset}', 'target_{onset}');
	xlabel('time (s)')
	ylabel('normalized position')
	nicegraph;
end

p = sqrt(HVF(:,1).^2+HVF(:,2).^2);
v = gradient(p,1/200);

sel = pl_confidence>0.8;
v(~sel) = NaN;
figure(43)
clf
ax(1) = subplot(311);
plot(tpl,p)
ax(2) = subplot(312);
plot(tpl,v)
linkaxes(ax,'x');


crit = 3*nanstd(v);
sel = v>crit;
selon                                       = [0;diff(sel)]; % onsets (+1)
onset = find(selon==1);
subplot(313)
plot(tpl,selon);
hold on
plot(tpl,sel);

% keyboard
subplot(312)
% idx = find(sel);
verline(tpl(onset));
horline(crit);

subplot(311)
verline(tpl(onset));
hold on
plot(tpl,pl_confidence);

%%
ntar = numel(tareve);

RT = NaN(ntar,1);
for ii = 1:ntar
	sel = tpl>tareve(ii) & tpl<(tareve(ii)+1);
% 	sum(sel)
	idx = find(sel);
	sel = onset>idx(1) & onset<idx(end);
	if sum(sel)
		rt = onset(sel)-idx(1);
		rt = rt/fs*1000; % ms
		RT(ii) = rt(1);
	end
end
sel = isnan(RT);
RT = RT(~sel);
figure(44)
clf
plotpost(RT);
%%


%% Resample to fixed 200 Hz
[HVF,t]				= resample(HVF,tpl,fs);
pl_confidence		= resample(pl_confidence,tpl,fs);

if dsp
	figure(1)
	clf
	plot(t,HVF,'linewidth',2);
	hold on
	
	plot(t,pl_confidence,'.','linewidth',2);
	ylim([-1.1 1.1]);
	verline(tareve,'k--');
	verline(fixeve,'g--');
	title('subject 002, session 005, 17-06-2022,');
	legend('horizontal','vertical','frontal','confidence','fixation_{onset}', 'target_{onset}');
	xlabel('time (s)')
	ylabel('normalized position')
	nicegraph;
end

%% Experimental block
idx			= repmat((-0.4*fs):(1*fs),ntrials,1); % which data per trial
idx			= round(idx+tareve'*fs); % index per time re target onset (note: events are not resampled)

% convert HVF to matrix NTRIALS x NSAMPLES
[H,V,F]		= deal(HVF(:,1),HVF(:,2),HVF(:,3));
H			= H(idx);
V			= V(idx);
F			= F(idx);

ti			= 1000*linspace(-0.4,1,length(H));
if dsp
	lH = lowpass(H','Fc',80,'Fs',fs,'order',20)';
	lV = lowpass(V','Fc',80,'Fs',fs,'order',20)';
	lF = lowpass(F','Fc',80,'Fs',fs,'order',20)';
	
	figure(2)
	clf
	
	subplot(231)
	plot(ti,H);
	hold on
	whos lH
	plot(ti,lH,'k-');
	subplot(232)
	plot(ti,V);
	
	subplot(233)
	plot(ti,F);
	
	
	subplot(234)
	plot(ti,gradient(lH,fs));
	
	subplot(235)
	plot(ti,gradient(lV,fs));
	
	subplot(236)
	plot(ti,gradient(lF,fs));
	
	str = {'horizontal','vertical','frontal'};
	for ii = 1:3
		subplot(2,3,ii)
		axis([min(ti) max(ti) -1 1]);
		nicegraph;
		xlabel('time (ms)');
		ylabel('normalized position');
		title(str{ii});
		
		subplot(2,3,ii+3)
		xlim([min(ti) max(ti)]);
		nicegraph;
		xlabel('time (ms)');
		ylabel('normalized velocity');
	end
end


% keyboard

%% Convert to dat and csv file
% create trial structure
[ntrials,nsamples]		= size(H);
nstim					= 3;

trial(ntrials).nstim	= nstim;

for ii = 1:ntrials
	trial(ii).nstim = nstim;
	for jj = 1:nstim
		switch jj
			case 1
				trial(ii).stim(jj).modality		= 'LED';
				trial(ii).stim(jj).azimuth		= 0;%  6) azimuth location
				trial(ii).stim(jj).elevation	= 0;%  6) azimuth location
				trial(ii).stim(jj).ondelay		= 0;%  6) azimuth location
				trial(ii).stim(jj).offdelay		= 300;%  6) azimuth location
				trial(ii).stim(jj).intensity	= 100;%  6) azimuth location
			case 2
				trial(ii).stim(jj).modality		= 'LED';
				trial(ii).stim(jj).azimuth		= tar(ii,1);%  6) azimuth location
				trial(ii).stim(jj).elevation	= tar(ii,2);%  6) azimuth location
				trial(ii).stim(jj).ondelay		= 300;%  6) azimuth location
				trial(ii).stim(jj).offdelay		= 1300;%  6) azimuth location
				trial(ii).stim(jj).intensity	= 100;%  6) azimuth location
			case 3
				trial(ii).stim(jj).modality		= 'data acquisition';
				trial(ii).stim(jj).azimuth		= [];%  6) azimuth location
				trial(ii).stim(jj).elevation	= [];%  6) azimuth location
				trial(ii).stim(jj).ondelay		= 0;%  6) azimuth location
				trial(ii).stim(jj).offdelay		= 300;%  6) azimuth location
				trial(ii).stim(jj).intensity	= [];%  6) azimuth location
		end
	end
end

% conversion
slims2hoop(H,V,F,fname,'001',trial);



%% Calibration block
idx			= repmat(-0.3*fs:0,ntrials,1);
ev = [fixeve(2:end) tareve(end)+1];
idx			= round(idx+ev'*fs);


[H,V,F]		= deal(HVF(:,1),HVF(:,2),HVF(:,3));
H			= H(idx);
V			= V(idx);
F			= F(idx);
ti			= 1000*linspace(-0.3,0,length(H));
if dsp
	lH = lowpass(H','Fc',80,'Fs',fs,'order',20)';
	lV = lowpass(V','Fc',80,'Fs',fs,'order',20)';
	lF = lowpass(F','Fc',80,'Fs',fs,'order',20)';
	
	figure(3)
	clf
	
	subplot(231)
	plot(ti,H);
	hold on
	whos lH
	plot(ti,lH,'k-');
	subplot(232)
	plot(ti,V);
	
	subplot(233)
	plot(ti,F);
	
	
	subplot(234)
	plot(ti,gradient(lH,fs));
	
	subplot(235)
	plot(ti,gradient(lV,fs));
	
	subplot(236)
	plot(ti,gradient(lF,fs));
	
	str = {'horizontal','vertical','frontal'};
	for ii = 1:3
		subplot(2,3,ii)
		axis([min(ti) max(ti) -1 1]);
		nicegraph;
		xlabel('time (ms)');
		ylabel('normalized position');
		title(str{ii});
		
		subplot(2,3,ii+3)
		xlim([min(ti) max(ti)]);
		nicegraph;
		xlabel('time (ms)');
		ylabel('normalized velocity');
	end
end


%% Convert to dat and csv file
% create trial structure
[ntrials,~]		= size(H);
nstim					= 2;

trial(ntrials).nstim	= nstim;

for ii = 1:ntrials
	trial(ii).nstim = nstim;
	for jj = 1:nstim
		switch jj
			case 1
				trial(ii).stim(jj).modality		= 'LED';
				trial(ii).stim(jj).azimuth		= tar(ii,1);%  6) azimuth location
				trial(ii).stim(jj).elevation	= tar(ii,2);%  6) azimuth location
				trial(ii).stim(jj).ondelay		= 0;%  6) azimuth location
				trial(ii).stim(jj).offdelay		= 300;%  6) azimuth location
				trial(ii).stim(jj).intensity	= 100;%  6) azimuth location
	
			case 2
				trial(ii).stim(jj).modality		= 'data acquisition';
				trial(ii).stim(jj).azimuth		= [];%  6) azimuth location
				trial(ii).stim(jj).elevation	= [];%  6) azimuth location
				trial(ii).stim(jj).ondelay		= 0;%  6) azimuth location
				trial(ii).stim(jj).offdelay		= 300;%  6) azimuth location
				trial(ii).stim(jj).intensity	= [];%  6) azimuth location
		end
	end
end

% conversion
slims2hoop(H,V,F,fname,'000',trial);

end

function slims2hoop(H,V,F,fname,block,trial)

%%
[~,name,~]	= fileparts(fname);
name			= [name '-' block];
ext             = '.dat';
datfile         = fcheckext(name,ext);
ext             = '.csv';
csvfile         = fcheckext(name,ext);
dlm				= ';';

%% Load data

[ntrials,nsamples]	= size(H);

%% Create dat-file
fid = fopen(datfile,'w','l');
for ii = 1:ntrials
	[x,y,z] = deal(H(ii,:)', V(ii,:)', F(ii,:)');
	d = [x y z zeros(length(x),5)]; % the configuration of the SPHERE coils is different from HOOP coils
	fwrite(fid,d,'float');
end
fclose(fid);


%% Create csv file


%% Header
Repeats			= 1;
ntrials_norep	= ntrials/Repeats;
Random			= 0;
nchan			= 8;
ITI				= [0 0];
Fs				= 200;
% first line, general configuration
first			= [0 ntrials_norep Repeats ntrials ITI Random nchan];
dlmwrite(csvfile,first,dlm)

% next lines, channel configuration
M				= NaN(nchan,6);
for ii			= 1:nchan
	M(ii,:)			= [0 ii ii round(Fs) round(Fs) nsamples];
end
dlmwrite(csvfile,M,'-append','delimiter',dlm)

%% Write trial information
fid = fopen(csvfile,'a+');
for trlIdx = 1:ntrials
	nstim = trial(trlIdx).nstim;
	for stmIdx = 1:nstim
		fprintf(fid,'%d%c',trlIdx,dlm); % 1) trial number
		fprintf(fid,'%d%c',stmIdx,dlm); % 2) stimulus number in trial
		fprintf(fid,'%d%c',round(mean(ITI)),dlm); % 3) Desired Inter Trial Interval
		fprintf(fid,'%d%c',round(mean(ITI)),dlm); % 4) Actual inter trial interval
		switch trial(trlIdx).stim(stmIdx).modality
			case 'LED'
				mod = 'led';
				attr = 0; % default
			case 'sound'
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
end




