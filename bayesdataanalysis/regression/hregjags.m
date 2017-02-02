function samples = hregjags(x,y,s)
% MCMC = HREGJAGS(X,Y,S)
%
% Example for Jags-Ymet-XmetSsubj-MrobustHier.R

%% Load data
if nargin<1
	% Load data file and specity column names of x (predictor) and y (predicted):
	file			= 'HierLinRegressData.csv';
	myData			= getcsvdata(file);
	xName			= 'X' ; yName = 'Y' ; sName='Subj';
	fileNameRoot	= 'HierLinRegressData-Jags-' ;
else
		xName			= 'X' ; yName = 'Y' ; sName='Subj';
	myData = struct(xName,x,yName,y,sName,s);
	fileNameRoot	= 'HierLinRegressData-Jags-' ;

end

% file='IncomeFamszState.csv';
% myData = csvread( file,1,0);
% xName			= 'Famsz' ; yName = 'Income' ; sName='State';
% fileNameRoot	= 'HIncomeFamszState-Lin-Jags-' ;

% file='BugsRatsData.csv';
% myData = csvread( file,1,0);
% xName = 'Day' ; yName = 'Weight' ; sName='Subj'
% fileNameRoot = 'BugsRatsData-Jags-'

graphFileType = 'eps';

%% Generate the MCMC chain:
%startTime = tic;
samples = genMCMC('data',myData,'xName',xName,'yName',yName,'sName',sName,...
	'numSavedSteps',2000,'thinSteps',1,'saveName',fileNameRoot);
%duration = toc(startTime);

% plotpost(samples.beta1mu)

%% Display diagnostics of chain, for specified parameters:
% for ( parName in c('beta0mu','beta1mu','nu','sigma','beta0[1]','beta1[1]') ) {
%  diagMCMC( codaObject=mcmcCoda , parName=parName ,
%            saveName=fileNameRoot , saveType=graphFileType )
% }
% %-------------------------------------------------------------------------------
% % Get summary statistics of chain:
% summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
% show(summaryInfo)
% % Display posterior information:
% plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , sName=sName ,
%           compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
%           pairsPlot=TRUE , showCurve=FALSE ,
%           saveName=fileNameRoot , saveType=graphFileType )
% %-------------------------------------------------------------------------------

%% Subfunctions
function samples = genMCMC(varargin)

chain_globals;
data			= keyval('data',varargin);
xName			= keyval('xName',varargin,'x');
yName			= keyval('yName',varargin,'y');
sName			= keyval('sName',varargin,'s');
numSavedSteps	= keyval('numSavedSteps',varargin,10000);
thinSteps		= keyval('thinSteps',varargin,1);
saveName		= keyval('saveName',varargin);
nChains			= keyval('nChains',varargin,nChainsDefault);

%% THE DATA.
y = data.(yName);
x = data.(xName);
% Convert sName to consecutive integers:
[~,~,s] = unique(data.(sName));
ns		= max(s);

% Do some checking that data make sense:
if any(~isfinite(y))
	errror('All y values must be finite.');
end
if any(~isfinite(x))
	error('All x values must be finite.');
end

%Ntotal = length(y)
% Specify the data in a list, for later shipment to JAGS:
dataStruct = struct('x',x,'y',y,'s',s,'Nsubj',ns);

writemodel;

%% RUN THE CHAINS
parameters = {'beta0','beta1','beta0mu','beta1mu','zbeta1sigma','zbeta0sigma','sigma','nu'};
burnInSteps		= 2000;			% Number of steps to 'burn-in' the samplers.
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethodDefault,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname = fullfile(pwd, 'model.txt');

initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).zbeta0			= zeros(ns,1);    % because data are standardized
	initsStruct(ii).zbeta1			= zeros(ns,1);        % because data are standardized
	initsStruct(ii).zbeta0mu		= 0;        % because data are standardized
	initsStruct(ii).zbeta1mu		= 0;        % because data are standardized
	initsStruct(ii).zbeta0sigma		= 1;        % because data are standardized
	initsStruct(ii).zbeta1sigma		= 1;        % because data are standardized
	initsStruct(ii).zsigma			= 1;        % because data are standardized
	initsStruct(ii).nuMinusOne		= 1;        % because data are standardized
end

% Use JAGS to Sample
samples = matjags( ...
	dataStruct, ...
	modelname, ...
	initsStruct, ...
	'doparallel' , doparallel, ...
	'nchains', nChains,...
	'nburnin', burnInSteps,...
	'nsamples', nIter, ...
	'thin', thinSteps, ...
	'monitorparams', parameters, ...
	'savejagsoutput' , 1 , ...
	'verbosity' , 1 , ...
	'dic',0,...
	'cleanup' , 0);

if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

samples.beta1sigma = samples.zbeta1sigma*std(y)/std(x);
samples.beta0sigma = samples.zbeta0sigma*std(y)+mean(y) + samples.zbeta0sigma*mean(x)*std(y)/std(x);


parameterNames = fieldnames(samples); % get all parameter names
for ii = 1:numel(parameterNames)
	n = size(samples.(parameterNames{ii}));
	
	switch numel(n)
		case 2
			samples.(parameterNames{ii}) = samples.(parameterNames{ii})(:);
		case 3
			samples.(parameterNames{ii}) = reshape(samples.(parameterNames{ii}),n(1)*n(2),n(3));
			
	end
end


% %===============================================================================
%
% smryMCMC = function(  codaSamples ,
%                       saveName=NULL ) {
%   mcmcMat = as.matrix(codaSamples,chains=FALSE)
%   paramNames = colnames(mcmcMat)
%   summaryInfo = NULL
%   for ( pName in paramNames ) {
%     summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
%   }
%   rownames(summaryInfo) = paramNames
%   if ( !is.null(saveName) ) {
%     write.csv( summaryInfo , file=paste(saveName,'SummaryInfo.csv',sep='') )
%   }
%   return( summaryInfo )
% }
%
% %===============================================================================
%
% plotMCMC = function( codaSamples , data , xName='x' , yName='y' , sName='s' ,
%                      compValBeta0=NULL , ropeBeta0=NULL ,
%                      compValBeta1=NULL , ropeBeta1=NULL ,
%                      compValSigma=NULL , ropeSigma=NULL ,
%                      showCurve=FALSE ,  pairsPlot=FALSE ,
%                      saveName=NULL , saveType='jpg' ) {
%   % showCurve is TRUE or FALSE and indicates whether the posterior should
%   %   be displayed as a histogram (by default) or by an approximate curve.
%   % pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
%   %   of parameters should be displayed.
%   %-----------------------------------------------------------------------------
%   y = data[,yName]
%   x = data[,xName]
%   s = factor(data[,sName])
%   nSubj = length(unique(s)) % should be same as max(s)
%   mcmcMat = as.matrix(codaSamples,chains=TRUE)
%   chainLength = NROW( mcmcMat )
%   beta0mu = mcmcMat[,'beta0mu']
%   beta1mu = mcmcMat[,'beta1mu']
%   sigma = mcmcMat[,'sigma']
%   zbeta0mu = mcmcMat[,'zbeta0mu']
%   zbeta1mu = mcmcMat[,'zbeta1mu']
%   zsigma = mcmcMat[,'zsigma']
%   nu = mcmcMat[,'nu']
%   log10nu = log10(nu)
%   %-----------------------------------------------------------------------------
%   if ( pairsPlot ) {
%     % Plot the parameters pairwise, to see correlations:
%     openGraph()
%     nPtToPlot = 1000
%     plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
%     panel.cor = function(x, y, digits=2, prefix='', cex.cor, ...) {
%       usr = par('usr'); on.exit(par(usr))
%       par(usr = c(0, 1, 0, 1))
%       r = (cor(x, y))
%       txt = format(c(r, 0.123456789), digits=digits)[1]
%       txt = paste(prefix, txt, sep='')
%       if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
%       text(0.5, 0.5, txt, cex=1.25 ) % was cex=cex.cor*r
%     }
%     pairs( cbind( beta0mu , beta1mu , sigma , log10nu )[plotIdx,] ,
%            labels=c( expression(mu[beta*0]) , expression(mu[beta*1]) ,
%                      expression(sigma) ,  expression(log10(nu)) ) ,
%            lower.panel=panel.cor , col='skyblue' )
%     if ( !is.null(saveName) ) {
%       saveGraph( file=paste(saveName,'PostPairs',sep=''), type=saveType)
%     }
%   }
%   %-----------------------------------------------------------------------------
%   % Marginal histograms:
%   % Set up window and layout:
%   nPtToPlot = 1000
%   plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
%   openGraph(width=8,height=8)
%   layout( matrix( 1:9 , nrow=3, byrow=TRUE ) )
%   par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
%   histInfo = plotPost( beta0mu , cex.lab = 1.75 , showCurve=showCurve ,
%                        compVal=compValBeta0 , ROPE=ropeBeta0 ,
%                        xlab=bquote(mu[beta*0]) , main=paste('Intercept, Group Level') )
%   histInfo = plotPost( beta1mu , cex.lab = 1.75 , showCurve=showCurve ,
%                        compVal=compValBeta1 , ROPE=ropeBeta1 ,
%                        xlab=bquote(mu[beta*1]) , main=paste('Slope, Group Level') )
%   plot( beta1mu[plotIdx] , beta0mu[plotIdx] ,
%         xlab=bquote(mu[beta*1]) , ylab=bquote(mu[beta*0]) ,
%         col='skyblue' , cex.lab = 1.75 )
%   histInfo = plotPost( zbeta0mu , cex.lab = 1.75 , showCurve=showCurve ,
%                        %compVal=compValBeta0 , ROPE=ropeBeta0 ,
%                        xlab=bquote(zmu[beta*0]) , main=paste('Intercept, Group Level') )
%   histInfo = plotPost( zbeta1mu , cex.lab = 1.75 , showCurve=showCurve ,
%                        %compVal=compValBeta1 , ROPE=ropeBeta1 ,
%                        xlab=bquote(zmu[beta*1]) , main=paste('Slope, Group Level') )
%   plot( zbeta1mu[plotIdx] ,zbeta0mu[plotIdx] ,
%         xlab=bquote(zmu[beta*1]) , ylab=bquote(zmu[beta*0]) ,
%         col='skyblue' , cex.lab = 1.75 )
%   histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
%                        compVal=compValSigma , ROPE=ropeSigma ,
%                        xlab=bquote(sigma) , main=paste('Scale, Subj Level') )
%   histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
%                        compVal=NULL , ROPE=NULL ,
%                        xlab=bquote(log10(nu)) , main=paste('Normality, Subj Level') )
%   plot( log10nu[plotIdx] , sigma[plotIdx] ,
%         xlab=bquote(log10(nu)) ,ylab=bquote(sigma) ,
%         col='skyblue' , cex.lab = 1.75 )
%   if ( !is.null(saveName) ) {
%     saveGraph( file=paste(saveName,'PostMarg',sep=''), type=saveType)
%   }
%   %-----------------------------------------------------------------------------
%   % Data with superimposed regression lines and noise distributions:
%   nPanels=25
%   nPlots = ceiling(nSubj/nPanels)
%   for ( plotIdx in 1:nPlots ) {
%     openGraph()
%     par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
%     layout(matrix(1:nPanels,nrow=5,byrow=TRUE))
%     xRang = max(x)-min(x)
%     yRang = max(y)-min(y)
%     xLimMult = 0.1
%     yLimMult = 0.1
%     xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
%     yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
%     %for ( sIdx in unique(ceiling(seq(1,nSubj,length=nPanels))) ) {
%     for ( sIdx in ((plotIdx-1)*nPanels+1):min(nSubj,(plotIdx-1)*nPanels+nPanels)) {
%       thisSrows = (as.numeric(s)==sIdx)
%       plot( x[thisSrows] , y[thisSrows] ,
%             cex=1.0 , lwd=1 , col='black' , xlim=xLim , ylim=yLim ,
%             xlab=xName , ylab=yName , cex.lab=1.0 ,
%             main=paste0('Unit: ',levels(s)[sIdx]) ,
%             cex.main=1.0  )
%       % Superimpose a smattering of believable regression lines:
%       nPredCurves=30
%       xComb = seq(xLim[1],xLim[2],length=301)
%       for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
%         b0 = mcmcMat[i,paste0('beta0[',sIdx,']')]
%         b1 = mcmcMat[i,paste0('beta1[',sIdx,']')]
%         lines( xComb , b0+b1*xComb , col='skyblue' )
%       }
%       points( x[thisSrows] , y[thisSrows] , pch=19 )
%     }
%     if ( !is.null(saveName) ) {
%       saveGraph( file=paste0(saveName,'PostPredSubj',plotIdx), type=saveType)
%     }
%   }
%   %-----------------------------------------------------------------------------
%   % Data with superimposed regression lines and noise distributions:
%   openGraph()
%   par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
%   % Plot data values:
%   xRang = max(x)-min(x)
%   yRang = max(y)-min(y)
%   xLimMult = 0.2
%   yLimMult = 0.2
%   xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
%   yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
%   plot( x , y , pch='' , cex=1.0 , col='black' ,
%         xlim=xLim , ylim=yLim ,
%         xlab=xName , ylab=yName , cex.lab=1.0 ,
%         main='All Units' , cex.main=1.0  )
%   % Superimpose a smattering of believable regression lines:
%   nPredCurves=70
%   for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
%     abline( mcmcMat[i,'beta0mu'] , mcmcMat[i,'beta1mu'] , col='skyblue' )
%   }
%   for ( sIdx in 1:nSubj ) {
%     thisSrows = (as.numeric(s)==sIdx)
%     lines( x[thisSrows] , y[thisSrows] , type='o' , pch=19 ) %, pch=sIdx , col=sIdx )
%   }
%   %
%   if ( !is.null(saveName) ) {
%     saveGraph( file=paste(saveName,'PostPredAll',sep=''), type=saveType)
%   }
%   %-----------------------------------------------------------------------------
%   % Individual parameter estimates (intercept, slope)
%   beta0hdi = apply( mcmcMat[, grep('^beta0\\[',colnames(mcmcMat)) ] , 2 , HDIofMCMC )
%   beta1hdi = apply( mcmcMat[, grep('^beta1\\[',colnames(mcmcMat)) ] , 2 , HDIofMCMC )
%   beta0range = range( beta0hdi )
%   beta1range = range( beta1hdi )
%   library('car') % install.packages('car')
%   openGraph()
%   par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
%   plot( 0 , 0 , type='n' ,
%         xlim=beta0range , ylim=beta1range ,
%         xlab='Intercept' , ylab='Slope' ,
%         main='Individual Parameters (mean and 95%CI)' )
%   for ( s in 1:nSubj ) {
%     b0s = mcmcMat[, grep(paste0('^beta0\\[',s,'\\]'),colnames(mcmcMat)) ]
%     b1s = mcmcMat[, grep(paste0('^beta1\\[',s,'\\]'),colnames(mcmcMat)) ]
%     dataEllipse(b0s,b1s,levels=c(0.95),plot.points=FALSE,col='skyblue',lwd=1)
%   }
%   if ( !is.null(saveName) ) {
%     saveGraph( file=paste(saveName,'PostSubjPar',sep=''), type=saveType)
%   }
% }
%
% %===============================================================================


function data = getcsvdata(file)
myData			= csvread( file,1,0);
fid		= fopen(file,'r');
str		= fgetl(fid);
fclose(fid);
idx		= strfind(str,',');
nIdx	= numel(idx)+1;
for ii = 1:nIdx
	switch ii
		case 1
			label = str(2:idx(ii)-2);
		case nIdx
			label = str(idx(ii-1)+2:end-1);
		otherwise
			label = str(idx(ii-1)+2:idx(ii)-2);
	end
	data.(label) = myData(:,ii);
end

function writemodel
str = ['data {\r\n',...
	'# Standardize the data:\r\n',...
	'\t\tNtotal <- length(y)\r\n',...
	'\t\txm <- mean(x)\r\n',...
	'\t\tym <- mean(y)\r\n',...
	'\t\txsd <- sd(x)\r\n',...
	'\t\tysd <- sd(y)\r\n',...
	'\t\tfor ( i in 1:length(y) ) {\r\n',...
	'\t\t\tzx[i] <- ( x[i] - xm ) / xsd\r\n',...
	'\t\t\tzy[i] <- ( y[i] - ym ) / ysd\r\n',...
	'\t\t}\r\n',...
	'}\r\n',...
	'# Specify the model for standardized data:\r\n',...
	'model {\r\n',...
	'\tfor ( i in 1:Ntotal ) {\r\n',...
	'\t\tzy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] , 1/zsigma^2 , nu )\r\n',...
	'\t}\r\n',...
	'\tfor ( j in 1:Nsubj ) {\r\n',...
	'\t\tzbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )\r\n',...
	'\t\tzbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )\r\n',...
	'\t}\r\n',...
	'# Priors vague on standardized scale:\r\n',...
	'zbeta0mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'zbeta1mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'zsigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'nu <- nuMinusOne+1\r\n',...
	'nuMinusOne ~ dexp(1/29.0)\r\n',...
	'# Transform to original scale:\r\n',...
	'for ( j in 1:Nsubj ) {\r\n',...
	'beta1[j] <- zbeta1[j] * ysd / xsd\r\n',...
	'beta0[j] <- zbeta0[j] * ysd  + ym - zbeta1[j] * xm * ysd / xsd\r\n',...
	'}\r\n',...
	'beta1mu <- zbeta1mu * ysd / xsd\r\n',...
	'beta0mu <- zbeta0mu * ysd  + ym - zbeta1mu * xm * ysd / xsd\r\n',...
	'sigma <- zsigma * ysd\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);
