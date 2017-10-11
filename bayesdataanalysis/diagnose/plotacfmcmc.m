function plotacfmcmc(mcmcStruct,varargin)
% PLOTACFMCMC(MCMC)
%
% Plot autocorrelation function of MCMC samples.
%
% See also DIAGMCMC, ESS

%% Initialization
parNames	= fieldnames(mcmcStruct);
parName		=  keyval('parName',varargin,parNames{1});

x			= transpose(mcmcStruct.(parName));
nChain		= size(x,2);
maxlag		= 35;

%% Autocorrelation
L = NaN(maxlag+1,nChain);
C = L;
for chainIdx = 1:nChain
	s				= x(:,chainIdx);
	s				= zscore(s);
	try % the Signal Processing Toolbox
		[acf,lag]		= xcorr(s,maxlag,'coeff');
		lag				= transpose(lag(maxlag+1:end));
		acf				= acf(maxlag+1:end);
	catch % or try another method
		lag		= -maxlag:maxlag;
		if ispc
		acf		= mxcorr(s,s,lag);
		else
			acf = ones(size(lag));
		end
		lag				= transpose(lag(maxlag+1:end));
		acf				= acf(maxlag+1:end);
	end
	L(:,chainIdx)	= lag;
	C(:,chainIdx)	= acf;
end

%% Effective sample size
e = ess(x);

%% Graphics
plot(L,C,'.-');
hold on
title(parName);
%   matplot( xMat , yMat , type='o' , pch=20 , col=plColors , ylim=c(0,1) ,
%            main='' , xlab='Lag' , ylab='Autocorrelation' )
%   abline(h=0,lty='dashed')
%   text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
%         labels=paste('ESS =',round(EffChnLngth,1)) )

xlim([-2 maxlag+2]);
ylim([-0.1 1.1]);
box off;
horline;
set(gca,'TickDir','out','XTick',0:5:35,'YTick',0:0.2:1);
xlabel('Lag');
ylabel('Autocorrelation');
str = ['ESS = ' num2str(e,'%0.0f')];
text(35,0.9,str,'HorizontalAlignment','right');
%   EffChnLngth = effectiveSize(mcmcStruct[,c(parName)])
