function bernouillitheta

close all;



%% Luuk van de Rijt's MATRIX Speech-In-Noise test Binomial
% cd('D:\SPIN-project\results_perSubjList\');
datadir = '/Volumes/mbaudit1/Marc van Wanrooij/SPIN/words';
cd(datadir);

groups		=  {'words','lists','subjects'};
ngroups		= numel(groups);
loadFlag	= true;
loadFlag = false;
sampleFlag	= [false false false];
% sampleFlag	= [true true true];


%% Load data
% First we load the data
if ~loadFlag
	tic
	d			= dir;
	d			= d([d.isdir]);
	dirnames	= {d(3:end).name}; % first 2 are . and ..;
	ndir		= numel(dirnames);
	T			= []; % presented words
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
		matfiles = dir('*.mat'); % find all mat-files
		matfiles = {matfiles.name};
		nfiles		= numel(matfiles);
		for fIdx = 1:nfiles % load every mat-file, containing data stored from experiment, = sentence
			fname = matfiles{fIdx}; 
			load(fname); % data from 1 sentence = 5 words
			rec

			stim		= rec.wordstimulus; % presented word
			res			= rec.wordresponse; % word response
			correct		= strcmp(stim,res); % determine correct response
			nstim		= numel(stim); % should be 5, as we presented 5 words
			trialNr		= rec.trialsCompleted; % number of trials/sentences completed per list, needed as every mat-file (1 sentence) contains speechlevels for all previous trials, and one in advance
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
			
			% check to see if correct responses as determined above
			% correspond to nCorrectTrial in rec structure
			sel = sum(correct)==rec.nCorrectTrial(trialNr);
			if ~sel
				disp('nCorrectTrial not corect');
				keyboard
			end
			
			
			if fIdx==2
				keyboard
			end
			
			snr			= repmat(speechlevel-noiselevel,1,nstim); % log(a)-log(b) = log(a/b)
			list		= repmat(rec.listNr,1,nstim); % list number (1-9)
			id			= rec.subjectId(1:end-2); % subject id
			subject		= repmat(str2double(id),1,nstim); % subject id times # words
			session		= repmat(str2double(rec.subjectId(end)),1,nstim); % subjects returned for up to 3 sessions
			modality	= repmat(mod,1,nstim);
			
			T			= [T stim]; %#ok<*AGROW> % words
			C			= [C correct];
			SNR			= [SNR snr];
			L			= [L list];
			S			= [S subject];
			O			= [O order];
			Z			= [Z session]; % zitting
			M			= [M modality];
		end
		cd ..
	end
	save spin T C SNR L S O Z M;
	toc
else
	load('spin');
end

keyboard
%% Sample data

[~,~,words] = unique(T); 	%% from strings in cell to numbers in double array
y			= double(C)';

sel			= M==2; % audio only, M==1
unique(SNR)
sel			= M==1 & SNR==-21; % audio only, M==1

y			= y(sel);

words		= words(sel)';
uwords		= unique(words);
nwords		= numel(uwords);

subjects	= S(sel);
usubjects	= unique(subjects);
nsubjects	= numel(usubjects);

mu			= NaN(nwords,nsubjects);
for jj = 1:nsubjects
	for ii = 1:nwords
		sel =	  words == uwords(ii) & subjects==usubjects(jj);
		mu(ii,jj) = mean(y(sel));
	end
end


%%
[muwords,idx]	= sort(median(mu,2));
musorted		= mu(idx,:); % by word performance



figure
subplot(121)
imagesc(mu(idx,:)')
axis square
set(gca,'YDir','normal');

subplot(122)
% plot(muwords,musorted,'k.');
plot(1:nwords,musorted,'k.');
hold on
axis square
% axis([-0.1 1.1 -0.1 1.1]);
% unityline;
box off

for ii = 1:nwords
	p = prctile(musorted(ii,:),[25 50 75]);
	P(ii,:) = p;
end

% X = muwords;
X = 1:nwords; X = X';
Y = smooth(P(:,2));
U = P(:,3)-Y;
L = Y-P(:,1);
% errorbar(X,Y,L,U,'ko-','MarkerFaceColor','w');
E = [smooth(P(:,1)),smooth(P(:,3))];
errorpatch(X',Y',E');

axis([-5 55 -0.1 1.1])
plot([-5 55],[-0.1 1.1],'r-','LineWidth',2);

%%
% keyboard
% hist(mu)
