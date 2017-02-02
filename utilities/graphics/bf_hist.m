function [h,f,xi] = bf_hist(x,varargin)
% H = BF_HIST(X)
%
% plot density of data X.
%
% Optional name-value input arguments include:
% - 'Color'
% - 'xlabel'
%
% see also KSDENSITY, HIST



%% Optional, undocumented:
col			= keyval('Color',varargin);
xlabels		= keyval('xlabel',varargin,'x');
lw			= keyval('LineWidth',varargin,2);
yoffset		= keyval('yoffset',varargin,0);
xi		= keyval('xi',varargin);
fun		= keyval('function',varargin,'pdf');

if isempty(xi)
[f,xi]		= ksdensity(x,'function',fun);
else
[f,xi]		= ksdensity(x,xi,'function',fun);
end
if ~isempty(col)
h			= plot(xi,f+yoffset,'-','Color',col,'LineWidth',lw);
else
	h			= plot(xi,f+yoffset,'-','LineWidth',lw); % let matlab decide
end
xlabel(xlabels);
axis square;
box off;
set(gca,'TickDir','out');

