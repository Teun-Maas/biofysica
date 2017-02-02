function plotloc(x,y,varargin)
% PLOTLOC(X,Y)
%
% Plot stimulus-response localization plot
%
% PLOTLOC(X,Y,NAME,VALUE,...)
% Additional NAME-VALUE pairs can be used
%	'az' - true (default) = azimuth text
%	'beta' - regression coefficients
%	'xlim'
%	'ylim'
%
% See also REGLINE, REGJAGS, REGSTATS

%% Initialization
azFlag		= keyval('az',varargin,'true');
xl			= keyval('xlim',varargin,[-100 100]);
yl			= keyval('ylim',varargin,[-100 100]);
xt			= keyval('xtick',varargin,-90:30:90);
yt			= keyval('ytick',varargin,-90:30:90);
col			= keyval('color',varargin,[.7 .7 .7]);
credFlag	= keyval('cred',varargin,true); % show credible regression lines
diagFlag	= keyval('diag',varargin,false);

beta		= keyval('beta',varargin);
if isempty(beta)
	if exist('matjags.m','file');
		beta = regjags(y,x,'diag',diagFlag);
	end
end

% r = corrcoef(x,y);

if isstruct(beta)
	b0		= beta.beta0(:);
	b1		= beta.beta1(:);
	mubeta = [mean(b0) mean(b1)];
	mur		= median(beta.pearson.r);
end

if diagFlag
	figure;
end
%% Plot data values:
hold on
xlim(xl);
ylim(yl);

%% Superimpose a smattering of believable regression lines:
if credFlag
	for ii =  round(linspace(1,length(b0),50))
		beta	= [b0(ii) b1(ii)];
		h		= regline(beta','k-');
		set(h,'Color',col);
	end
end
h = regline(mubeta','k-');
set(h,'LineWidth',2);
plot(x,y,'ko','MarkerFaceColor',col,'LineWidth',1,'MarkerSize',8);

%% Graphics
axis square;
xlim(xl);
ylim(yl);
set(gca,'TickDir','out','XTick',xt,'YTick',yt);
box off
xlabel('Target (deg)');
ylabel('Response (deg)');
title('Azimuth Data with credible regression lines');
unityline('k:');
s	= sign(round(mubeta(1)))+2;
sgn = {'-','+','+'};
if azFlag
	greek = 'alpha';
else
	greek = 'epsilon';
end
str = {['\' greek '_R = ' num2str(mubeta(2),'%.1f')  '\' greek '_T ' sgn{s} ' ' num2str(abs(mubeta(1)),'%.0f')],...
	['R^2=' num2str(mur.^2,'%.2f')]};

h = text(xl(2)*0.95,yl(1)*0.8,str,'HorizontalAlignment','right','FontSize',12,'FontWeight','bold','Color','r');
