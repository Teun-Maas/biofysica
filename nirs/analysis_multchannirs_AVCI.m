function analysis_multchannirs_AVCI(varargin)

close all;
warning off;
%% Initialization
gr = 'CI';
d				= keyval('dir',varargin,['/Users/marcw/DATA/Roos Cartignij/NIRS sessie' filesep gr]);
cd(d);

%% Do stuff
expdir			= dir;
expnames		= {expdir.name};
idxO2				= vectorstrfind(expnames,'.');
expnames(idxO2)	= [];
nexp			= numel(expnames);
close all

% for expIdx = [1 3:nexp]
% for expIdx = 1:nexp
for expIdx = 5
	
	expname = expnames{expIdx};
	cd(d);
	cd(expname)
	fname	= dir('rcs*.mat');
	fname	= fname.name;
	load(fname)
	
	keyboard
	[data_epoched,data_avg] = epoch(data_rcs,eventstim,nevents,uevent,uevents,fname);
	title(fname);
	
	drawnow
	
	%% GLM
	trial	= data_rcs.trial{1};
	ntrial	= size(trial,1);
	X		= data_rcs.HDR;
	whos X
% 	X = X(1:11,:);
	T = [];
	P = [];
	for ii = 1:ntrial
		Y = trial(ii,:);
		b = regstats(Y',X','linear',{'tstat'});
		T = [T b.tstat.t];
		P = [P b.tstat.pval];
	end
	
	load('ciavlay.mat')
	figure(901)
	clf
	figure(902)
	clf
	ColMap = statcolor(64,[],[],[],'def',8);
	for ii = 2:(nevents+1)
		row = uevent(ii-1,1)+1;
		col = uevent(ii-1,2)+1;
		sb = (row-1)*4+col;
		
		label		= data_rcs.label;
		triallabel	= data_avg(ii-1).label;
		
		idx			= vectorstrfind(label,'[O2Hb]');
		label		= label(idx);
		trial		= data_avg(ii-1).avg;
		
		trial		= trial(idx,:);
		[seldat, sellay] = match_str(label,lay.label);
		
		pos			= lay.pos(sellay,:);
		tstat		= T(ii,seldat)';
		trial		= trial(seldat,:);
		pval		= P(ii,seldat)';
		
		for jj = 1:2
			switch jj
				case 1
					sel		= pos(:,1)<0;
				case 2
					sel		= pos(:,1)>0;
			end
			X		= pos(sel,1);
			Y		= pos(sel,2);
			L		= label(sel);
			for lidx = 1:numel(L)
				L{lidx} = L{lidx}(1:end-7);
			end
			Tr		= trial(sel,:);
			Z		= tstat(sel);
			Zp		= pval(sel);
			F		= scatteredInterpolant(X,Y,Z,'natural');
			Fp		= scatteredInterpolant(X,Y,Zp,'natural');
			Xi		= linspace(min(X),max(X),100);
			Yi		= linspace(min(Y),max(Y),100);
			[XI,YI] = meshgrid(Xi,Yi);
			Zi		= F(XI,YI);
			Zpi		= Fp(XI,YI);
			% 	Zi(Zpi>0.05) = NaN;
			figure(900+jj)
			colormap(ColMap);
			
			subplot(4,4,sb)
			
			imagesc(Xi,Yi,Zi);
			hold on
% 			plot(X,Y,'ws','MarkerSize',10);
			
			%%
			for kk = 1:size(Tr,1)
				y = Tr(kk,:);
				y = y-y(50);
				x = X(kk)+(1:length(y))-length(y)/2;
				plot(x,repmat(Y(kk),size(x)),'w-')
				plot([x(50) x(50)],[min(500*y+Y(kk)) max(500*y+Y(kk))],'w-')

				plot(x,500*y+Y(kk),'k-')
			end
			%%
			set(gca,'YDir','normal',...
				'YTick',[],'XTick',[]);
			axis square
			axis([min(X)-50 max(X)+50 min(Y)-50 max(Y)+50]);
			box off
			caxis([-15 15]);
			if (ii-1)<12
				title(ueventstr(ii-1,:))
			end
			colorbar
		end
	end
	cd('/Users/marcw/DATA/Roos Cartignij/matlab');
	
	drawnow
	pause(.5)
	
% 	%% Epoching
% 	% 	close all
% 	trg			= data_rcs.trigger{1}; % trigger indices for onset of the standard
% 	ntrg		= numel(trg);
% 	eventstim		= eventstim(1:ntrg);
% 	for ii = 1:nevents
% 		cfg						= [];
% 		cfg.dataset				= fname;
% 		cfg.trialdef.eventvalue = 1; % read conditions
% 		cfg.trialdef.prestim    = 5; % in seconds
% 		cfg.trialdef.poststim   = 10; % in seconds
% 		
% 		
% 		sel						= eventstim==uevents(ii); % auditory 75%?
% 		
% 		trig					= trg(sel);
% 		sevents					= eventstim(sel);
% 		trl						= NaN(length(trig)-2,4);
% 		
% 		for stdidx = 1:length(trig)
% 			begsample		= trig(stdidx) - cfg.trialdef.prestim*data_rcs.fsample;
% 			endsample		= trig(stdidx) + cfg.trialdef.poststim*data_rcs.fsample - 1;
% 			offset			= -cfg.trialdef.prestim*data_rcs.fsample;
% 			trigger			= sevents(stdidx); % remember the trigger (=condition) for each trial
% 			trl(stdidx, :)	= [round([begsample endsample offset])  trigger];
% 		end
% 		
% 		cfg.trl			= trl;
% 		cfg.lpfilter	= 'yes';                              % apply lowpass filter
% 		cfg.lpfreq		= 0.3;                                 % lowpass
% 		cfg.hpfilter   = 'yes';                              % apply highpass filter
% 		cfg.hpfreq     = 0.1;                              %
% 		data_tst = data_rcs;
% 		data_epoched = ft_redefinetrial(cfg,data_tst);
% 		
% 		
% 		cfg				= [];
% 		
% 		data_avg		= ft_timelockanalysis(cfg, data_epoched);
% 		% 	data_avg		= ft_timelockanalysis(cfg, data_avg);
% 		
% 		
% 		load('ciavlay.mat')
% 		
% 		figure(100)
% 		cfg				= [];
% 		cfg.showlabels	= 'yes';
% 		cfg.layout		= lay;
% 		cfg.interactive = 'yes';
% 		cfg.channel		= '* [O2Hb]';
% 		% 		cfg.channel		= '* [HHb]';
% 		
% 		cfg.ylim		= [-0.2 0.2];
% 		cfg.baseline	= [-5 0];
% 		
% 		%cfg.graphcolor = 'r';
% 		ft_multiplotER(cfg, data_avg)
% 		pause
% 	end
	
end
%%
% keyboard

function [data_epoched,data_avg] = epoch(data,events,nevents,uevent,uevents,fname)
%% Epoching
close all
trg			= data.trigger{1}; % trigger indices for onset of the standard
ntrg		= numel(trg);
events		= events(1:ntrg);
for ii = 1:nevents
	cfg						= [];
	cfg.dataset				= fname;
	cfg.trialdef.eventvalue = 1; % read conditions
	cfg.trialdef.prestim    = 5; % in seconds
	cfg.trialdef.poststim   = 10; % in seconds
	
	
	sel						= events==uevents(ii); % auditory 75%?
	trig					= trg(sel);
	sevents					= events(sel);
	trl						= NaN(length(trig)-2,4);
	
	for stdidx = 1:length(trig)
		begsample		= trig(stdidx) - cfg.trialdef.prestim*data.fsample;
		endsample		= trig(stdidx) + cfg.trialdef.poststim*data.fsample - 1;
		offset			= -cfg.trialdef.prestim*data.fsample;
		trigger			= sevents(stdidx); % remember the trigger (=condition) for each trial
		trl(stdidx, :)	= [round([begsample endsample offset])  trigger];
	end
	
	cfg.trl			= trl;
	% 	cfg.lpfilter	= 'yes';                              % apply lowpass filter
	% 	cfg.lpfreq		= 0.3;                                 % lowpass
	% 	cfg.hpfilter   = 'yes';                              % apply highpass filter
	% 	cfg.hpfreq     = 0.1;
	cfg.baseline	= [-5 0];
	%
	data_epoched	= ft_redefinetrial(cfg,data);
	
	cfg				= [];
	
	data_avg(ii)		= ft_timelockanalysis(cfg, data_epoched);
	
	label = data_avg(ii).label;
	trial = data_avg(ii).avg;
	idx		= vectorstrfind(label,'[O2Hb]');
	trial = trial(idx,:);
	time = data_avg(ii).time;
	ntrial = size(trial,1);
	
end

