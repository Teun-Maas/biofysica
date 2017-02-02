function psifitexample

%% Simulation
datFlag = 4;
[x,y,s] = getdata(datFlag);

% keyboard
close all;
cd /Users/marcw/Documents/MATLAB

switch datFlag
	case 2
		centroid = 'median';
		guessrate = 'infer';
		lapserate = 'infer';
		fun = @logisticfun;
	case 6
		centroid = 'mean';
		guessrate = 0.5;
		lapserate = 'infer';
		fun = @revweibullfun;
	case 7
		centroid = 'mode';
		guessrate = 0;
		lapserate = 0;
		fun = @logisticfun;
	case 4
				centroid = 'mean';
		guessrate = 0.1;
		lapserate = 0;
		fun = @logisticfun;
	otherwise
		centroid = 'median';
		guessrate = 'infer';
		lapserate = 'infer';
		fun = @logisticfun;
end

t = tic;
psifit(x,y,s,'gamma',guessrate,'lambda',lapserate,'function',fun,'numSavedSteps',25000,'thinSteps',3,'showCentroid',centroid,'showDiag',true);
e = toc(t)

savegraph(mfilename,'png');
function [x,y,s] = getdata(datFlag)


switch datFlag
	case 1
		%% DBDA data
		fname	= 'HtWtData110.csv';
		fname	= which(fname);
		myData	= csvread(fname,1,0);
		y		= myData(:,1);
		x		= myData(:,3);
		s		= ones(size(x));
	case 2
		%% Bayesian Causal Modeling workshop Michael Lee & Eric
		% load data
		p = which('data_x.txt'); % does it exist?
		[pathstr,~,~] = fileparts(p);
		fname	= 'data_x.txt';
		x		= dlmread(fullfile(pathstr,fname),'\t');
		fname	= 'data_n.txt';
		n		= dlmread(fullfile(pathstr,fname),'\t');
		fname	= 'data_r.txt';
		r		= dlmread(fullfile(pathstr,fname),'\t');
		nsubjs	= size(x,1);
		% 'transform' data
		S		= repmat(transpose(1:nsubjs),1,size(x,2)); % subject vector
		X		= x(:);
		Y		= [r(:) n(:)]; % 2-dimensional y containing rate and number of trials
		S		= S(:);
		% remove NaNs
		sel		= ~isnan(x);
		x		= X(sel);
		y		= Y(sel,:);
		s		= S(sel,:);
		
	case 3
		%% Luuk van de Rijt's MATRIX Speech-In-Noise test Binomial
		cd('/Users/marcw/DATA/Test/SPIN');
		d = dir('List*');
		nlists = numel(d);
		X = [];
		Y = [];
		S = [];
		L = [];
		for lstIdx = 1:nlists
			cd(d( lstIdx).name)
			fnames = dir('*.mat');
			nfiles = numel(fnames);
			x = NaN(nfiles,20);
			y = x;
			s = x;
			r = x;
			l = x;
			for fIdx = 1:nfiles
				fname = fnames(fIdx).name;
				load(fname);
				x(fIdx,:) = rec.speechLevels(1:20);
				y(fIdx,:) = rec.nCorrectTrial(1:20);
				s(fIdx,:) = str2num(fname(1));
				r(fIdx,:) = str2num(fname(3));
				l(fIdx,:) = str2num(fname(end-4));
			end
			
			X = [X; x(:)];
			S = [S; s(:)];
			Y = [Y; y(:)];
			L = [L; l(:)];
			cd ..
		end
		
		x = X;
		y = Y;
		s = S;
		l = L;
		
		n = numel(x);
		
		Y = [];
		X = [];
		S = [];
		L = [];
		for ii = 1:n
			X = [X; repmat(x(ii),5,1)]; %#ok<*AGROW>
			Y = [Y; ones(y(ii),1); zeros(5-y(ii),1)];
			S = [S; repmat(s(ii),5,1)]; %#ok<*AGROW>
			L = [L; repmat(l(ii),5,1)]; %#ok<*AGROW>
			
		end
		x = X;
		y = Y;
		s = S;
		s = L;
	case 4
		%% Luuk van de Rijt's MATRIX Speech-In-Noise test Bernouilli
		dr = '/Users/marcw/Gitlab/biofysica/bayesdataanalysis/data/spin';
		cd(dr);
		d			= dir;
		dirnames	= {d(4:end).name}; % first 3 are . and .. and .DS_Store;
		ndir		= numel(dirnames);
		T			= [];
		L			= [];
		S			= [];
		C			= [];
		O			= [];
		SNR			= [];
		exc			= {'koopt','negen','vuile','boeken',''}; % error in code
		order		= 1:5;
		
		for dirIdx = 1:ndir
			% 				dirnames{dirIdx};
			cd(dirnames{dirIdx});
			matfiles = dir('*.mat');
			matfiles = {matfiles.name};
			nfiles = numel(matfiles);
			for fIdx = 1:nfiles
				fname = matfiles{fIdx};
				load(fname);
				stim		= rec.wordstimulus;
				if strcmp(stim,exc)
					% 						matfiles{fIdx}
					stim = {'Jan','koopt','negen','vuile','boeken'};
				end
				res			= rec.wordresponse;
				correct		= strcmp(stim,res);
				nstim		= numel(stim);
				
				% bug - speechlevel is not recorded per trial
				% try to find in fname
				tr		= fname;
				idx1	= strfind(tr,'_tr');
				idx2	= strfind(tr,'_');
				if idx2(2)>idx1
					idx2	= idx2(2)-1;
				else
					idx2	= idx2(3)-1;
				end
				idx1	= idx1+3;
				tr		= str2double(tr(idx1:idx2));
				level	= rec.speechLevels(tr);
				
				
				snr			= repmat(level-str2double(rec.settings.noiseLvl),1,nstim);
				list		= repmat(rec.listNr,1,nstim);
				subject		= repmat(str2double(rec.subjectId(1)),1,nstim);
				T			= [T stim];
				C			= [C correct];
				SNR			= [SNR snr];
				L			= [L list];
				S			= [S subject];
				O			= [O order];
			end
			cd ..
		end
		
		[~,~,t] = unique(T);
		
		x = SNR';
		y = double(C)';
		% 		s = S'; % subjects
		% 		s = L'; % lists
		
		s	= t; % words
		sel = ismember(s,1:5);
		x	= x(sel);
		y	= y(sel);
		s	= s(sel);
		% 		s = O'; % word order in sentence
	case 5
		% Sebastian Ausili's ITD/ILD JND
		fname = which('ild02_SA.mat');
		load(fname)
		x = uml.x;
		y = uml.r;
		
		% add some X=0 Y=random Bernouilli values
		n = 50;
		x = [zeros(n,1);x];
		y = [rndval(0,1,[n,1]);y];
		s = ones(size(x));
		
	case 6
		fname = which('ILD_const_high_300.mat');
		load(fname)
		
		
	case 7
		fname = which('ILD_const_low_ABS_300d.mat');
		
		load(fname)
		y = yr;
		
		[x,~,subs] = unique(x);
		r = accumarray(subs,y,[],@sum);
		n				= accumarray(subs,ones(size(subs)),[],@sum);
		y = [r n];
		s = ones(size(x));
		
end

