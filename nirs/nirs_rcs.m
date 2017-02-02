function data_rcs = nirs_rcs(data,varargin)
% DATA = NIRS_RCS(DATA)
%
% Perform reference channel subtraction on DATA, by regressing out all
% shallow/reference channels from each deep channel.
%
% A shallow channel is defined as having a distance between transmitter and
% receiver less than 3 cm, deep channels have a T-R distance >= 3cm.
%
% To do:
% - remove motion/initialization/gain setting artefacts (now first
% 4000 samples are removed).
% - Check whether multiple linear regression makes sense (i.e. do far
% reference channels have less influence?)
% - check whether MLR with both HbO and HbR reference channels makes sense
%
% See also NIRS_OPTODELAYOUT, NIRS_OPTODEDISTANCE

%% Initialization

if nargin<1
	close all
	w			= what('LR-04-2015-06-17_active');
	DataFolder	= w.path;
	cd(DataFolder)
	fname		= 'data_trans.mat';
	load(fname);
	data		= data_trans;
end

dispFlag	= keyval('disp',varargin,false);
method		= keyval('method',varargin,'closest');

%% Relevant parameters
label	= data.label; % transformed channel label
% flabel	= data.opto.fiberlabel; % fiber label
% label	= label(1:2:end); % deoxy labeling removed

% %% Correcting a flaw in nirs
% for tIdx = 1:16
% 	flabel{end-16+tIdx} = ['Tx' num2str(tIdx)]; % number from 1 to 16
% end

%% Obtaining receiver-transmitter fiber distances
if ~isfield(data,'chanpos')
	data = nirs_optodelayout(data);
end
if ~isfield(data.opto,'fiberdistance')
	data	= nirs_optodedistance(data);
end
d			= data.opto.fiberdistance;
pos			= data.chanpos;

%% Distance and position
dat			= data.trial{:};
% dat			= dat(1:2:end,:); % oxy

selshallow	= d~=max(d);
seldeep		= d==max(d);

shallow		= dat(selshallow,:);
% shallow		= bsxfun(@minus,shallow,mean(shallow,2)); % mean detrend
shallowpos	= pos(selshallow,:);
shallowlabel = label(selshallow);

deep		= dat(seldeep,:);
% deep		= bsxfun(@minus,deep,mean(deep,2)); % mean detrend
deeppos		= pos(seldeep,:);
deeplabel	= label(seldeep);

time		= data.time{:};

%% Reference channel subtraction
ndeep		= size(deep,1);
signal		= NaN(size(deep));
x			= shallow';
sb			= ceil(sqrt(ndeep));
for dpIdx	= 1:ndeep
	y				= deep(dpIdx,:)';
	l				= deeplabel{dpIdx}(end-5:end);
	sell			= vectorstrfind(shallowlabel,l,'logic',true);
	
	if strcmpi(method,'closest');
		p				= deeppos(dpIdx,:); % current deep position
		d				= sqrt( (p(1)-shallowpos(:,1)).^2+(p(2)-shallowpos(:,2)).^2 ); % current distance from shallow channels
		mnd				= min(d); % find closest shallow/reference channels
		seld			= d==mnd;
		
		
		sel				= sell & seld;
	else
		sel = sell;
	end
	b				= regstats(y,x(:,sel),'linear',{'beta','r'});
	res				= b.r;
	signal(dpIdx,:) = res';
	
	if dispFlag
		subplot(sb,sb,dpIdx)
		plot(y)
		hold on
		plot(res);
		hold off
		xlim([0 length(y)]);
		str = [data.label{dpIdx}];
		bf_text(0.1,0.9,str);
		title(str)
		axis off
		
		drawnow
	end
	
end

%% Auto-correlation
% rshallow	= corrcoef(shallow'); % correlation matrix shallow channels
% rdeep		= corrcoef(deep'); % correlation matrix deep channels
% rsignal		= corrcoef(signal'); % correlation matrix signal channels

%% overwrite
data_rcs			= data;
data_rcs.time{1}	= time;
data_rcs.trial{1}	= signal;
data_rcs.label		= label(seldeep);
data_rcs.chanpos	= data_rcs.chanpos(seldeep,:);
data_rcs.sci		= data_rcs.sci(seldeep);

% %% Some more graphics
% figure;
% subplot(131)
% imagesc(rshallow);
% % colorbar;
% caxis([0 1]);
% set(gca,'YDir','normal','XTick',1:sum(selshallow),'YTick',1:sum(selshallow),'TickDir','out',...
% 	'XTickLabel',label(selshallow),'XTickLabelRotation',90,...
% 	'YTickLabel',label(selshallow));
% axis square;
% title('Shallow')
%
% subplot(132)
% imagesc(rdeep);
% % colorbar;
% caxis([0 1]);
% set(gca,'YDir','normal','XTick',1:sum(seldeep),'YTick',1:sum(seldeep),'TickDir','out',...
% 	'XTickLabel',label(seldeep),'XTickLabelRotation',90,...
% 	'YTickLabel',label(seldeep));
% axis square;
% title('Deep')
%
% subplot(133)
% imagesc(rsignal);
% % colorbar;
% caxis([0 1]);
% set(gca,'YDir','normal','XTick',1:sum(seldeep),'YTick',1:sum(seldeep),'TickDir','out',...
% 	'XTickLabel',label(seldeep),'XTickLabelRotation',90,...
% 	'YTickLabel',label(seldeep));
% axis square;
% title('Signal')
%
% figure
% plot(rdeep.^2,rsignal.^2,'k.');
% axis square;
% ylim([0 1]);
% pa_unityline;
% box off
% set(gca,'TickDir','out','Xtick',0:0.2:1,'YTick',0:0.2:1);
% xlabel('R^2 Deep')
% ylabel('R^2 Signal');
%
% %% Even more graphics
% figure;
% subplot(311)
% plot(time,shallow');
% ylim([-15 15]);
% title('Shallow channels')
% xlabel('time (s)');
% ylabel('\DeltaHbO_2');
% box off
% set(gca,'TickDir','out');
%
% subplot(312)
% plot(time,deep');
% ylim([-15 15]);
% title('Deep channels')
% xlabel('time (s)');
% ylabel('\DeltaHbO_2');
% box off
% set(gca,'TickDir','out');
%
% subplot(313)
% plot(time,signal');
% ylim([-15 15]);
% title('Signal channels')
% xlabel('time (s)');
% ylabel('\DeltaHbO_2');
% box off
% set(gca,'TickDir','out');


% pa_datadir;
% print('-dpng','-r300',mfilename);

% function rcs(dat)
