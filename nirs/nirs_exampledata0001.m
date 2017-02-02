function nirs_exampledata0001

%% Initialization

close all
clear variables

%% Load mat files
cd('/Users/marcw/Gitlab/biofysica/nirs/example_data/exampledata0001'); %#ok<*UNRCH> % contains all relevant data files

d		= dir('nirs0*.mat');
[m,~]	= size(d);
nmbr	= NaN(m,1);
for ii	= 1:m
	fname		= d(ii).name;
	nmbr(ii)	= str2double(fname(5:7));
	
end
[~,Txt]			= xlsread('exampledata0001.xlsx');
Txt				= Txt(2:end,:);
figure(2)
clf
subj	= 17;
sel		= nmbr == subj;
fnames	= {d(sel).name};
nfiles	= numel(fnames);
chan	= 2;
S		= [];
So		= [];
T		= [];
M		= [];
E		= [];
for jj = 1:nfiles
	load(fnames{jj});
	S	= [S;nirs.signal]; %#ok<*NODEF,*AGROW>
	So	= [So;nirs.deepchan]; %#ok<*AGROW>
	%% add timings, continuous
	if ~isempty(E)
		e = T(end)*nirs.Fs;
		t = T(end);
		T = [T;nirs.time+t];
		E = [E [nirs.event.sample]+e];
	else
		T = [T;nirs.time];
		E = [E [nirs.event.sample]];
	end
	fname	= fnames{jj}
	fname	= fname(5:end-4)
	sel		= strcmpi(fname,Txt(:,1));
	
	%% Check for interleaved blocks
	modal	= Txt(sel,5);
	
	if strncmp(modal,'Random',5);
		M = [M {nirs.event.stim}];
	elseif strncmp(modal,'Audiovisual',11);
		a = {nirs.event.stim};
		for kk = 1:size(a,2)
			a {kk} = 'AV';
		end
		M = [M a];
	elseif strncmp(modal,'Auditory',8);
		a = {nirs.event.stim};
		for kk = 1:size(a,2)
			a {kk} = 'A';
		end
		M = [M a];
	elseif strncmp(modal,'Visual',6);
		a = {nirs.event.stim};
		for kk = 1:size(a,2)
			a {kk} = 'V';
		end
		M = [M a];
	end
end
S = pa_zscore(S);

%% Block average
clear N
N.event.sample	= E';
N.event.stim	= M;
N.Fs			= nirs.Fs;
N.fsdown		= nirs.fsdown;
col = {'r';'g';'b'};

mod			= {'A';'AV';'V'};
mx = NaN(3,1);
for chanIdx = chan
	for modIdx	= 1:3
		MU	= pa_nirs_blockavg(N,S(:,chanIdx),mod{modIdx});
		x			= 1:length(MU);
		x			= x/10;
		baseline	= nanmean(MU(:,100),2);
		MU2			= bsxfun(@minus,MU,baseline);
		
		figure(2)
		subplot(122)
		plot(x,nanmean(MU2),'-','Color',col{modIdx},'LineWidth',2)
		hold on
		
		pa_errorpatch(x,nanmean(MU2),nanstd(MU)./sqrt(size(MU,1)),col{modIdx});
		mx(modIdx) = max(nanmean(MU2));
	end
end
axis square;
box off
title(ii);

%% Hemo
nirs	= N;
fs		= nirs.Fs;
fd		= nirs.fsdown;
R		= S(:,1);
on		= ceil([nirs.event.sample]*fd/fs);
off		= on(2:2:end);
on		= on(1:2:end);
stim	= nirs.event.stim;
ustim	= unique(stim);
nstim	= numel(ustim);
n		= numel(R);
HDR		= NaN(3,n);
for sIdx = 1:nstim
	sel		= strcmp(ustim(sIdx),stim(1:2:end)) | strcmp('AV',stim(1:2:end));
	ons		= on(sel);
	offs	= off(sel);
	N		= length(R);
	bl      = zeros(1, N);
	for kk	= 1:length(ons)
		bl(ons(kk):offs(kk)) = 1;
	end
	hemo	= pa_nirs_hdrfunction(1,bl);
	hemo	= removeheartbeat([hemo' hemo'],0.1); % heartbeat
	hemo	= removeheartbeat(hemo,0.1,2); % respiration (0.2 Hz)
	hemo	= removeheartbeat(hemo,0.1,[0.05 .2]); % Mayer Wave[0.05 .2]
	hemo	= hemo(:,1);
	flim	= [0.008 .1];
	hemo	= pa_bandpass(hemo,flim,5); % we bandpass between 0.016 and 0.8 Hz to remove noise and cardiac component
	hemo	= hemo./max(hemo);
	HDR(sIdx,:) = hemo;
end


%% GlM
on		= ceil([nirs.event.sample]*fd/fs);
off		= on(2:2:end);
on		= on(1:2:end);
stim	= nirs.event.stim;
stim	= stim(1:2:end);
ustim	= unique(stim);
% disp('===================')
for chanIdx = chan
	y		= S(:,chanIdx);
	
	x		= 1:length(y);
	x		= x/10;
	on		= ceil([nirs.event.sample]*fd/fs);
	off		= on(2:2:end);
	on		= on(1:2:end);
	subplot(121)
	for kk = 1:numel(on)
		sel = x>on(kk)/10 & x<off(kk)/10+20;
		y1 = y(sel);
		y1 = y1-y1(1);
		plot(x(sel)-on(1)/10,y1,'k-');
		hold on;
	end
	
x		= HDR';
	y		= S(:,chanIdx);
	[b,~,stats]		= glmfit(x,y,'normal'); % A AV V

	b
	stats.t>0
	stats.t
	stats.p
	
	b(1)	= 0;
	figure(2)
	yfit	= glmval(b,x,'identity');
	yfit(isnan(y)) = NaN;
	x		= 1:length(yfit);
	x		= x/10;
	x = x-on(1)/10;
	sel = x>0 & x<650;
	subplot(121)
	plot(x(sel),yfit(sel),'Color',[.7 .7 .7])
	hold on

	box off
	axis square
	

end
% keyboard
subplot(121)
ax	= axis;
for ii = 1:numel(on)
	idx = strcmp(stim(ii),ustim);
	subplot(121)
	hp	= patch(([on(ii) on(ii) off(ii) off(ii)]-on(1))/10,[ax(3) ax(4) ax(4) ax(3)],col{idx});
	set(hp,'EdgeColor','none');
	alpha(hp,0.3);
end

%%
subplot(121)
xlim([-100 800])
ylim([-1 3]);
set(gca,'XTick',0:200:1000,'YTick',0:2,'TickDir','out',...
	'TickLength',[0.005 0.025],...
	'FontSize',10,'FontAngle','Italic','FontName','Helvetica'); % publication-quality
xlabel('time re stimulus #1 (s)','FontSize',12);
ylabel('\Delta HbO_2 (au)','FontSize',12);
title([]);
axis square;
box off
h = pa_text(0.05,0.9,'A');
set(h,'FontSize',15,'FontWeight','bold');
str = sprintf('Subject %d, Right channel, block 1 (A, V, AV interleaved)',subj);
title(str,'FontSize',12);
pa_horline(0,'k-')
% legend('Trace','Best-fit GLM','Location','SE');
% pa_horline(mx,'k:')

subplot(122)
ylim([-1 3]);
xlim([-5 50])
pa_errorpatch([10 30],[-6 -6],[20 20]);
set(gca,'XTick',0:10:40,'XTickLabel',(0:10:40)-10,'YTick',0:2,'TickDir','out',...
	'TickLength',[0.005 0.025],...
	'FontSize',10,'FontAngle','Italic','FontName','Helvetica'); % publication-quality
xlabel('time re stimulus onset (s)','FontSize',12);
ylabel('Average \Delta HbO_2 (au)','FontSize',12);
title([]);
axis square;
box off
h = pa_text(0.05,0.9,'B');
set(h,'FontSize',15,'FontWeight','bold');
str = sprintf('Subject %d, Right channel, block 1 (A, V, AV interleaved)',subj);
title(str,'FontSize',12);
pa_horline(0,'k-')
% pa_horline(mx,'k:')

legend('Visual','Audiovisual','Auditory','Location','SE');

%%
pa_datadir;
print('-depsc','-painters',[mfilename 'HbO2']);

