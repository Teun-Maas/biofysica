function plotdensmcmc(mcmcStruct,varargin)
% PLOTDENSMCMC(MCMC)
%
% Plot density of MCMC sample chains.
%
% See also DIAGMCMC

%% Initialization
parNames	= fieldnames(mcmcStruct);
parName		=  keyval('parName',varargin,parNames{1});

x			= transpose(mcmcStruct.(parName));
nChain		= nchain(x);

xMat = NaN(nChain,100);
yMat = xMat;
hdiLims = NaN(nChain,2);
for  cIdx = 1:nChain
	[yMat(cIdx,:),xMat(cIdx,:)] = ksdensity(x(:,cIdx));
	hdiLims(cIdx,:) = hdimcmc(x(:,cIdx));
end

plot(xMat',yMat','-');
hold on

box off;
set(gca,'TickDir','out');
xlabel('Parameter value');
ylabel('Density');
horline;

verline( hdiLims,'k:');
text( mean(hdiLims(:)), 0,'95% HDI','HorizontalAlignment','center','VerticalAlignment','bottom');
EffChnLngth = ess(x);
MCSE = std(x(:))/sqrt(EffChnLngth);
text( max(xMat(:)),max(yMat(:)),['MCSE = ' num2str(MCSE,'%.3f')],'HorizontalAlignment','center','VerticalAlignment','bottom');

