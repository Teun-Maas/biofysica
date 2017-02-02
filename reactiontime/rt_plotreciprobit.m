function rt_plotreciprobit(fname,varargin)
% RT_PLOTRECIPROBIT

% 2015 Marc van Wanrooij
% e: m.vanwanrooij@donders.ru.nl

% close all
% clear all
% clc

printFlag	= keyval('print',varargin,false);
col			= keyval('Color',varargin,'r');
xl			= keyval('xlim',varargin,[]);
xt			= keyval('xtick',varargin,[]);

if isstring(fname)
	load(fname);
elseif isnumeric(fname)
	RT = fname;
elseif isstruct(fname)
	RT = fname.RT;
end

if ~any(RT>100)
	RT = RT*1000;
	% Reaction times in ms
end

% pa_datadir;
% load('reactiontime')
% load('reactiontime.mat'); % load the reaction time data into Matlab workspace


%% Probit
% raw data

x = -1./sort((RT)); % multiply by -1 to mirror abscissa
n = numel(RT); % number of data points
y = probit((1:n)./n); % cumulative probability for every data point converted to prbt scale
y(isinf(y)) = probit(0.99);
plot(x,y,'o','Color','k','MarkerFaceColor',col);
hold on

% quantiles
p		= [1 2 5 10:10:90 95 98 99]/100;
prbt	= probit(p);
hold on

if isempty(xt)
xtick	= sort(-1./(150+[0 oct2bw(50,-1:5)])); % some arbitrary xticks
else
	xtick = xt;
end
set(gca,'XTick',xtick,'XTickLabel',-1./xtick);

if isempty(xl)
xlim([min(xtick) max(xtick)]);
else
	xl = [-1/xl(1) -1./xl(2)];
xlim(xl);
end	
set(gca,'YTick',prbt,'YTickLabel',p*100);
ylim([probit(0.1/100) probit(99.9/100)]);
axis square;
box off
xlabel('Reaction time (ms)');
ylabel('Cumulative probability (%)');

% this should be a straight line

sel = isfinite(y);
% b = regstats(y(sel),x(sel),'linear','beta');
% b = b.beta;
% h = regline(b,'k-');
% set(h,'Color',col,'LineWidth',1);
% 
b = robustfit(x(sel),y(sel));
h = regline(b,'k--');
set(h,'Color',col,'LineWidth',1);

horline(0,'k:');
if printFlag
	% optional, save figure
	savegraph(mfilename,'eps');
end
