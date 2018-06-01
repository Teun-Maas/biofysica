function [trial,cfg] = readexp(cfg)
% [TRIAL,CFG] = READEXP(FNAME)
%
% Load TRIAL and ConFiGuration parameters from exp-file FNAME
%
% See also GETEXP

% 2015 Marc van Wanrooij

%% Initialization
if nargin<1
	fname = fcheckexist([]);
	cfg = [];
else
	if isstruct(cfg)
		fname = fullfile(cfg.expdir,cfg.expfname);
	else
		fname	= cfg;
		cfg		= [];
	end
end

%% Open
fid		= fopen(fname,'r');

%% Skip comment lines
% the header of the exp file can include any number of comment lines.
% Typically it is best to give a clear description of the experiment in
% these lines
comment			= checkcomment(fid);
cfg.comment		= comment;

%% Header lines
% The header lines of the exp file contain some parameters on the overall
% experiment, such as # trials, # repeats, etc.
% Most are obsolete:
% - ITI / Inter trial interval: should typically be zero, as the experiment
%		should not take unnecessary time
% - Trials / # trials: obtained from actual number of trials
% - Repeats / # repeats: if you want to have multiple repetitions, include that in
%			the trial section of the exp file (it is better to have control
%			over the experiment, than to instruct a algorithm to take
%			control)
% - Random / randomization: take control and randomize the trials with
%			RANDPERM
% - Motor / power on of motor: the SPHERE setup does not have a motor
%			installed (yet)
% All header lines can be followed by a comment (which are ignored).
%
% Sphere exp files are allowed to have any number of header lines, but any
% lines not referring to the above parameters should abide by the
% convention:
% NAME INTEGER
% NAME should not include any spaces
% NAME and INTEGER should be separated by space or tab
% INTEGER should be an integer....
cfg			= getheader(fid,cfg);
% warnings
if isfield(cfg,'Lab')
	if isfield(cfg,'ITI')
		if any(cfg.ITI~=0)
			message = 'Do you really want to include ''ITI'' ?';
			warning(message);
		end
	end
	if isfield(cfg,'Repeats')
		message = 'Variable ''Repeats'' is not used';
		warning(message);
	end
	if isfield(cfg,'Random')
		message = 'Variable ''Random'' is not used';
		warning(message);
	end
	if isfield(cfg,'Motor')
		message = 'Variable ''Motor'' is not used in SPHERE lab';
		warning(message);
	end
else
	cfg.Lab = 1; % Hoop = default
	% 2 = Sphere
	% 3 = SphereMinor
	% 4 = NIRS-EEG?
end

%% Trial

%			MOD	X	Y	ID	INT	On		On		Off		Off		Event;	| 10
%						edg	bit	Event	Time	Event	Time
% LED		+	+	+	+	+	+		+		+		+				| 7
% SKY		+	+	+		+	+		+		+		+				| 7
% SKY		+	+		+	+	+		+		+		+				| 7
% SND		+	+	+	+	+	+		+								| 6
% ACQ		+					+		+								| 2
% LAS		+	+	+	+	+	+		+		+		+				| 6
% TRG0		+			+	+	+		+						+		| 5
% INP1/2	+															| 1

trial				= struct([]);
% to order fields, create them first:
trial(1).stim(1).X			= [];
trial(1).stim(1).Y			= [];
trial(1).stim(1).channel	= [];
trial(1).stim(1).detect		= [];
trial(1).stim(1).event		= [];
trial(1).stim(1).intensity	= [];
trial(1).stim(1).modality	= [];
trial(1).stim(1).offdelay	= [];
trial(1).stim(1).offevent	= [];
trial(1).stim(1).ondelay	= [];
trial(1).stim(1).onevent	= [];

tn		= 0;
isLas	= false; % if there are lasers in the exp-file, give a warning only once
while ~feof(fid)
	curLine             = fgetl(fid);
	firstCell           = sscanf(curLine,'%s',1);
	nchar				= length(firstCell);
	switch upper(firstCell)
		case '==>'
			tn	= tn+1; % trial number count
			sn	= 0; % stimulus number count
		case {'SND1','SND2','SND'} % could also have 7
			sn								= sn+1;
			par		= sscanf(curLine(nchar+1:end),'%d%d%d%f%d%d',[6,1]);
			trial(tn).stim(sn).modality		= 'sound';
			trial(tn).stim(sn).X			= par(1);
			trial(tn).stim(sn).Y			= par(2);
			trial(tn).stim(sn).matfile		= ['snd' num2str(par(3),'%03i') '.mat'];
			trial(tn).stim(sn).wavfile		= ['snd' num2str(par(3),'%03i') '.wav'];	 % double
			trial(tn).stim(sn).intensity	= par(4);
			trial(tn).stim(sn).onevent		= par(5);
			trial(tn).stim(sn).ondelay		= par(6);
			trial(tn).stim(sn).offevent		= par(5); % default duration
			trial(tn).stim(sn).offdelay		= par(6)+150; % default duration
		case 'LED'
			if	ismember(cfg.Lab,[1 4]) % Is this correct for sphereMinor?
				sn								= sn+1;
				par	= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d%d',[7,1]);
				trial(tn).stim(sn).modality		= 'LED';
				trial(tn).stim(sn).X			= par(1);
				trial(tn).stim(sn).Y			= par(2);
				trial(tn).stim(sn).intensity	= par(3);  % hoop: range 0-255, sphere range 1-50
				trial(tn).stim(sn).onevent		= par(4);
				trial(tn).stim(sn).ondelay		= par(5);
				trial(tn).stim(sn).offevent		= par(6);
				trial(tn).stim(sn).offdelay		= par(7);
			elseif	ismember(cfg.Lab,[2 3]) % sphere and sphereMinor with LED colours
				sn								= sn+1;
				par	= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d%d',[8,1]);
				trial(tn).stim(sn).modality		= 'LED';
				trial(tn).stim(sn).X			= par(1);
				trial(tn).stim(sn).Y			= par(2);
				trial(tn).stim(sn).colour		= par(3); % 0 - red, 1 - green
				trial(tn).stim(sn).intensity	= par(4);  % hoop: range 0-255, sphere range 1-50
				trial(tn).stim(sn).onevent		= par(5);
				trial(tn).stim(sn).ondelay		= par(6);
				trial(tn).stim(sn).offevent		= par(7);
				trial(tn).stim(sn).offdelay		= par(8);
			end
		case 'IRLED'
			sn								= sn+1;
			par	= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d%d%d',[8,1]);
			trial(tn).stim(sn).modality		= 'LED';
			trial(tn).stim(sn).X			= par(1);
			trial(tn).stim(sn).Y			= par(2);
			trial(tn).stim(sn).color		= par(3);
			trial(tn).stim(sn).intensity	= par(4);
			trial(tn).stim(sn).onevent		= par(5);
			trial(tn).stim(sn).ondelay		= par(6);
			trial(tn).stim(sn).offevent		= par(7);
			trial(tn).stim(sn).offdelay		= par(8);
		case 'SKY'
			sn								= sn+1;
			if ~sscanf(curLine(nchar+1:end),'%d',1) % center SKY led
				par	= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d%d',[7,1]);
				trial(tn).stim(sn).color			= par(2);
				trial(tn).stim(sn).Y			= 0; % bug: hoopxy2azel needs input
			else % not center sky
				par	= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d%d%d',[8,1]);
				trial(tn).stim(sn).Y			= par(2);
			end
			trial(tn).stim(sn).modality		= 'sky';
			trial(tn).stim(sn).X			= par(1);
			trial(tn).stim(sn).intensity	= par(3);
			trial(tn).stim(sn).onevent		= par(4);
			trial(tn).stim(sn).ondelay		= par(5);
			trial(tn).stim(sn).offevent		= par(6);
			trial(tn).stim(sn).offdelay		= par(7);
		case 'TRG0'
			sn								= sn+1;
			par			= sscanf(curLine(nchar+1:end),'%d%d%d%d%d',[5,1]);
			trial(tn).stim(sn).modality		= 'trigger';
			if par(1) == 1
				trial(tn).stim(sn).detect		= 'rise';
			elseif par==2
				trial(tn).stim(sn).detect		= 'fall';
			end
			trial(tn).stim(sn).channel		= par(2);
			trial(tn).stim(sn).onevent		= par(3);
			trial(tn).stim(sn).ondelay		= par(4);
			trial(tn).stim(sn).event		= par(5);
		case 'LAS'
			% does not exist yet for sphere, assume fixation LED
			sn = sn+1;
			par		= sscanf(curLine(nchar+1:end),'%d%d%d%d%d%d',[6,1]);
			if ~isLas
				disp('LAS detected; assuming a fixation LED')
				isLas = true;
			end
			
			trial(tn).stim(sn).modality		= 'LED';
			trial(tn).stim(sn).channel		= par(1);
			trial(tn).stim(sn).intensity	= par(6);
			trial(tn).stim(sn).onevent		= par(2);
			trial(tn).stim(sn).ondelay		= par(3);
			trial(tn).stim(sn).offevent		= par(4);
			trial(tn).stim(sn).offdelay		= par(5);
			trial(tn).stim(sn).X			= 0;
			trial(tn).stim(sn).Y			= 12;
		case 'ACQ'
			sn = sn+1;
			par		= sscanf(curLine(nchar+1:end),'%d%d',[2,1]); % could also be 3
			trial(tn).stim(sn).modality		= 'data acquisition';
			trial(tn).stim(sn).onevent		= par(1);
			trial(tn).stim(sn).ondelay		= par(2);
		case {'INP1','INP2','INP'}
			sn = sn+1;
			trial(tn).stim(sn).modality		= 'sound acquisition';
	end
end


%% check number of trials
ntrials			= numel(trial);
if isfield(cfg,'Trials') % Trials is the default name for nuber of trials from HOOP setup
	cfg.ntrials		= cfg.Trials; % we want cfg.ntrials
end
if isfield(cfg,'ntrials') % check whether number of trials in exp header corresponds to actual # of trials in exp file
	sel				= cfg.ntrials == ntrials;
	if ~sel
		message = 'Number of actual trials does not match number of trials stated in header';
		warning(message);
		disp(['CFG.NTRIALS set to ' num2str(ntrials)]);
		cfg.ntrials = numel(trial);
	end
else % or if it does not exist, set cfg.ntrials to # trials
	cfg.ntrials = ntrials;
	disp(['CFG.NTRIALS set to ' num2str(ntrials)]);
end

%% Extra: Azimuth and Elevation from Hoop
if cfg.Lab==2
	cfg		= spherelookup(cfg);
elseif cfg.Lab==4
	cfg		= spherelookupMinor(cfg);
end
for trlIdx = 1:cfg.ntrials % for every trial
	s			= trial(trlIdx).stim;
	trial(trlIdx).nstim = numel(s); % number of stimuli per trial
	for stmIdx	= 1:trial(trlIdx).nstim % for every stimulus in a trial
		X			= trial(trlIdx).stim(stmIdx).X;
		Y			= trial(trlIdx).stim(stmIdx).Y;
		mod			= trial(trlIdx).stim(stmIdx).modality;
		if ~isempty(X) % for every stimulus that has an X and Y parameter, determine azimuth and elevation
			if cfg.Lab==1 % Hoop lab
				if strcmpi(mod,'sky')
					[Az,El] = hoopsky2azel(X,Y);
				else
					[Az,El]	= hoopXY2azel(X,Y);
				end
			elseif ismember(cfg.Lab,[2 3]) % Sphere lab
				channel = cfg.interpolant(X,Y);
				Az		= cfg.lookup(channel+1,5);
				El		= cfg.lookup(channel+1,6);
			elseif ismember(cfg.Lab,4) % SphereMinor lab
				channel = cfg.interpolant(X,Y);
				Az		= cfg.lookup(channel+1,4);
				El		= cfg.lookup(channel+1,5);
			end
			trial(trlIdx).stim(stmIdx).azimuth		= Az;
			trial(trlIdx).stim(stmIdx).elevation	= El;
		end
	end
end

%% Close
fclose(fid);

function comment = checkcomment(fid)
isComment	= true;
cnt			= 0; % counter
comment		= cell(1);
while isComment % do this for every line that starts with '%'
	position	= ftell(fid); % find the position in the file
	str			= fscanf(fid,'%s',1); % read the string (moving the position in the file)
	commentline = fgetl(fid); % get the entire line (again repositioning)
	isComment	= strncmp(str,'%',1); % and check whether the first string of the line actually indicated a comment
	if ~isempty(commentline) && isComment
		cnt				= cnt+1;
		comment(cnt)	= {commentline};
	end
end
fseek(fid,position,'bof'); % let's get back to the position that was not a comment

function cfg = getheader(fid,cfg)

% isTrial			= false;
% while ~isTrial
% 	position	= ftell(fid);
% 	str			= fscanf(fid,'%s',1);
% 	isTrial	= strcmp(str,'==>');
% end
% fseek(fid,position,'bof');

cnt			= 0;
isTrial		= false;
header		= cell(1);
while ~isTrial
	position	= ftell(fid);
	str			= fscanf(fid,'%s',1);
	isTrial		= strcmp(str,'==>');
	if ~isTrial
		fseek(fid,position,'bof');
		cnt			= cnt+1;
		header(cnt)	= {fscanf(fid,'%s',1)};
		switch lower(header{cnt})
			case 'iti' % Inter-trial interval
				cfg.(header{cnt})	= fscanf(fid,'%d %d',[2 1]); % 2 integers: minimum and maximum possible inter trial interval
			case 'motor' % Motor
				cfg.(header{cnt})	= fscanf(fid,'%s',1); % String: yes or no
			otherwise
				cfg.(header{cnt})	= fscanf(fid,'%d',1); % Integer
		end
		checkcomment(fid);
	end
end
fseek(fid,position,'bof');

cfg.header		= header;
