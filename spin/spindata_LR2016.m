function spindata_LR2016


%% Luuk van de Rijt's MATRIX Speech-In-Noise test Binomial
% cd('D:\SPIN-project\results_perSubjList\');
datadir = '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);



%% Load data
% First we load the data
tic
d			= dir;
d			= d([d.isdir]);
dirnames	= {d(3:end).name}; % first 2 are . and ..;

ndir		= numel(dirnames);
T			= []; % presented words
R			= []; % response words
L			= []; % presented list
S			= []; % subject
C			= []; % correct response (0 or 1)
O			= []; % presentatio  order
SNR			= []; % Signal-to-noise ratio (speech-noise)
Z			= []; % session number (Z = Dutch zitting)
M			= []; % Modality of word
order		= 1:5;

for dirIdx = 1:ndir % loop through every data-directory = subject
	cd(dirnames{dirIdx});
	matfiles	= dir('*.mat'); % find all mat-files
	matfiles	= {matfiles.name};
	
	% remove trial 0: turning on experiment also creates a trial-mat
	idx			= vectorstrfind(matfiles,'tr0');
	sel			= ~ismember(1:numel(matfiles),idx);
	matfiles	= matfiles(sel);
	
	%
	nfiles		= numel(matfiles);
	for fIdx = 1:nfiles % load every mat-file, containing data stored from experiment, = sentence
		fname = matfiles{fIdx};
		load(fname); % data from 1 sentence = 5 words
		stim		= rec.wordstimulus; % presented word
		res			= rec.wordresponse; % word response
		sentence	= rec.trialsCompleted+20*(rec.listNr-1);
		correct		= strcmp(stim,res); % determine correct response
		nstim		= numel(stim); % should be 5, as we presented 5 words
		trialNr		= rec.trialsCompleted; % number of trials/sentences completed per list, needed as every mat-file (1 sentence) contains speechlevels for all previous trials
		speechlevel	= rec.speechLevels(trialNr); % this is changed in experiment for Luuk vd Rijt 2016
		noiselevel	= rec.noiseLevels(trialNr); % this is constant for Luuk vd Rijt 2016
		mod			= rec.settings.testtype; % Modality 'Audio+Video', 'Audio only', 'Video only'
		switch mod % change text to number
			case 'Audio + video'
				mod = 3;
			case 'Audio only'
				mod = 1;
			case 'Video only'
				mod = 2;
		end
		
		
		snr			= repmat(speechlevel-noiselevel,1,nstim); % log(a)-log(b) = log(a/b)
		list		= repmat(rec.listNr,1,nstim); % list number (1-9)
		id			= rec.subjectId(1:end-2); % subject id
		sentence	= repmat(sentence,1,nstim);
		% correct error for subject 2
		if str2double(id)==99
			id = fname(1);
		end
		
		%
		subject		= repmat(str2double(id),1,nstim); % subject id times # words
		session		= repmat(str2double(rec.subjectId(end)),1,nstim); % subjects returned for up to 3 sessions
		modality	= repmat(mod,1,nstim);
		
		T			= [T stim]; %#ok<*AGROW> % words
		R			= [R res];
		C			= [C correct];
		SNR			= [SNR snr];
		L			= [L list];
		Z			= [Z sentence];
		S			= [S subject];
		O			= [O order];
% 		Z			= [Z session]; % zitting
		M			= [M modality];
	end
	cd ..
end
save spin T R C SNR L S O Z M;
toc

