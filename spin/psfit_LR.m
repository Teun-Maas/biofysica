function psfit_LR

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
if ~loadFlag
	tic
	d			= dir;
	d = d([d.isdir]);
	dirnames	= {d(3:end).name}; % first 2 are . and ..;
	ndir		= numel(dirnames);
	T			= [];
	L			= [];
	S			= [];
	C			= [];
	O			= [];
	SNR			= [];
	Z			= [];
	M			= [];
	order		= 1:5;
	
	for dirIdx = 1:ndir
		cd(dirnames{dirIdx});
		matfiles = dir('*.mat');
		matfiles = {matfiles.name};
		nfiles = numel(matfiles);
		for fIdx = 1:nfiles
			fname = matfiles{fIdx};
			load(fname);
			stim		= rec.wordstimulus;
			res			= rec.wordresponse;
			correct		= strcmp(stim,res);
			nstim		= numel(stim);
			trialNr		= rec.trialsCompleted;
			speechlevel		= rec.speechLevels(trialNr);
			noiselevel		= rec.noiseLevels(trialNr);
			mod = rec.settings.testtype;
			switch mod
				case 'Audio + video'
					mod = 3;
				case 'Audio only'
					mod = 1;
				case 'Video only'
					mod = 2;
			end
			
			
			snr			= repmat(speechlevel-noiselevel,1,nstim); % log(a)-log(b) = log(a/b)
			list		= repmat(rec.listNr,1,nstim);
			id			= rec.subjectId(1:end-2);
			subject		= repmat(str2double(id),1,nstim);
			session		= repmat(str2double(rec.subjectId(end)),1,nstim);
			modality		= repmat(mod,1,nstim);
			
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
for grpIdx = 1:ngroups
	sel			= strcmp(groups,groups{grpIdx});
	sample		= sampleFlag(sel);
	group		= groups{sel};
	switch group
		case 'words'
			[~,~,t] = unique(T); 	%% from strings in cell to numbers in double array
		case 'lists'
			[~,~,t] = unique(L);
		case 'subjects'
			[~,~,t] = unique(S);
	end
	
	%% AV
	x		= SNR';
	y		= double(C)';
	s		= t;
	sel		= M==3; % audio only, M==1
	xAV		= x(sel);
	yAV		= y(sel);
	sAV		= s(sel);
	
	if sample
		samplesAV = psifit(xAV,yAV,sAV,'showPred','false'); %#ok<*USENS>
	end
	
	%% A
	x		= SNR';
	y		= double(C)';
	s		= t;
	sel		= M==1; % audio only, M==1
	xA		= x(sel);
	yA		= y(sel);
	sA		= s(sel);
	
	if sample
		samplesA = psifit(xA,yA,sA,'showPred','false'); %#ok<*UNRCH>
		save(['spinsample' group],'samplesA','samplesAV')
	else
		load(['spinsample' group])
	end
	
	%%
	close all
	figure(200)
	posteriorprediction(xAV,yAV,sAV,samplesAV,@logisticfun,'mean',[0 0.7 0]); %#ok<*NODEF>
	posteriorprediction(xA,yA,sA,samplesA,@logisticfun,'median',[0.6 0.6 1]);
	
	
	%% 
	parameterNames	= fieldnames(samplesA); % get all parameter names
	sel				= strfind(parameterNames,'mu');
	sel				= ~cellfun(@isempty,sel);
	parameterNames	= parameterNames(sel);
	npar			= numel(parameterNames);
	muA				= NaN(npar,1);
	for parIdx		= 1:npar
		muA(parIdx) = mean(samplesA.(parameterNames{parIdx}));
	end
	xWid		= max(xA)-min(xA);
	xComb		= linspace(min(xA)-0.1*xWid,max(xA)+0.1*xWid,201);
	ypred		= psifun(xComb,muA(1),muA(2),muA(3),muA(4),0.1,'function',@logisticfun);
	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
	
	%% V
	us = unique(s);
ns = length(us);
sb = ceil(sqrt(ns));

	y		= double(C)';
	s		= t;
	sel		= M==2; % audio only, M==1
	yV		= y(sel);
	sV		= s(sel);
	usV		= unique(sV);
	nsV		= numel(usV);
	mu		= NaN(nsV,1);
	for ii = 1:nsV
		sel = usV(ii) == sV;
		mu(ii) = mean(yV(sel));
		
		subplot(sb,sb,ii)
		horline(mu(ii),'r-');
		set(gca,'XTick',unique(xA),'XTickLabel',[],'YTick',0:0.2:1,'YTickLabel',[])
		
		pred = mu(ii)+ypred-mu(ii)*ypred;
		plot(xComb,pred,'k:','LineWidth',2)
	end
	
	%% Samples
	figure(300)
	clf
	parameterNames	= fieldnames(samplesA); % get all parameter names
	sel = strfind(parameterNames,'mu');
	sel = ~cellfun(@isempty,sel);
	parameterNames = parameterNames(sel);
	for parIdx			= 1:numel(parameterNames)
		a =[samplesA.(parameterNames{parIdx}); samplesAV.(parameterNames{parIdx})]';
		xax = minmax(a);
		figure(300)
		subplot(1,4,parIdx)
		plotpost(samplesA.(parameterNames{parIdx}),'showCurve',true,'xlim',xax,'xlab',parameterNames{parIdx}(3:end),'col',[0.5 0.5 1]);
		plotpost(samplesAV.(parameterNames{parIdx}),'showCurve',true,'xlim',xax,'xlab',parameterNames{parIdx}(3:end),'col',[0 0.7 0]);
		axis square
	end
	
	subplot(1,4,3);
	y		= double(C)';
	sel		= M==2; % audio only, M==1
	yV		= y(sel);
	verline(mean(yV),'r-');
	
	%%
	drawnow
	
	figure(200)
	print('-depsc',['individual' group mfilename])
	
	figure(300)
	print('-depsc',['group' group mfilename])
end
%%
function posteriorprediction(x,y,s,samples,fun,centroidFlag,col)
% POSTERIORPREDICTION(X,Y,B)
%
%

%% Parameters
theta			= samples.theta;
omega			= samples.omega;
gamma			= samples.gamma;
lambda			= samples.lambda;

% guess		= samples.guess;

%% Determine centroid
if strcmp(centroidFlag,'mode')
	centroid = ['bf_' centroidFlag]; % for random samples, BF_MODE works better than MODE
else
	centroid = centroidFlag;
end
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	str = [parameterNames{parIdx} ' = samples.(parameterNames{parIdx});'];
	eval(str);
	str = ['c' parameterNames{parIdx} ' = ' centroid '(samples.(parameterNames{parIdx}));'];
	eval(str);
end

%%
chainLength = length(theta);

%% believable psychometric curve
cVec		= floor(linspace(1,chainLength,30));
xWid		= max(x)-min(x);
xComb		= linspace(min(x)-0.1*xWid,max(x)+0.1*xWid,201);
xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}));
	selw		= xComb>=0;
	xComb		= xComb(selw);
end

us = unique(s);
ns = length(us);
sb = ceil(sqrt(ns));
sbcnt	= 0;


for ii = 1:ns % for every subject/group
	sbcnt		= sbcnt+1; % new subplot count
	sel			= s==us(ii);
	
	subplot(sb,sb,sbcnt)
	hold on
	
	%% Posterior credible psychometric curves
	for cIdx	= cVec
		ypred	= psifun(xComb,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',col);
	end
	
	
	%% Predictive posterior max?
	ypred = psifun(xComb,ctheta(ii),comega(ii),cgamma(ii),clambda(ii),0.1,'function',fun);
	
	
	%% Graphics
	plot(xComb,ypred,'k-','LineWidth',2,'Color',0.5*col);
	
	% Data
	if size(y,2)==1 % y = bernouilli rate
		[ux,~,subs]		= unique(x(sel));
		r				= accumarray(subs,y(sel),[],@sum);
		n				= accumarray(subs,ones(size(subs)),[],@sum);
	elseif size(y,2)==2 % y = [rate n]
		r = y(sel,1);
		n = y(sel,2);
		ux = x(sel);
	end
	rprop			= r./n;
	
	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5);
	
	
	%% Labels
	set(gca,'TickDir','out','XTick',unique(x),'YTick',0:0.2:1);
	axis square;
	box off;
	ylim([-0.1 1.2]);
	xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
	bf_text(0.1,0.9,num2str(us(ii)));
end
