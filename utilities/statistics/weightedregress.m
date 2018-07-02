function [betai,sei,XI,Ni] = weightedregress(X,Y,Z,sigma,XI,varargin)
% [BETA,SE,XI] = WEIGHTEDREGRESS(X,Y,Z,SIGMA,XI)
%
% Perform linear regression  Y vs Z weighted by factor X.
% Weighing is performed either through a Gaussian or boxcar profile with a
% width of SIGMA, centered  around XI.
%
% [BETA,SE,XI] = WEIGHTEDREGRESS(...,'NBIN',NBIN)
% Optionally, you can provide the minimum number NBIN of data points required
% within plusmin 1 SIGMA of XI (default = 0). This is useful when few
% data-points ara available at XI, and regression might be biased towards
% outliers.
%
% [BETA,SE,XI] = WEIGHTEDREGRESS(...,'WFUN',WFUN)
% Optionally, you can provide weight function: 'gaussian' or 'boxcar'
% (default WFUN = 'gaussian').
%
% [BETA,SE,XI] = WEIGHTEDREGRESS(...,'WFUN',WFUN)
% Optionally, you can provide the fit function: 'regstats' or 'robustfit'
%
% [BETA,SE,XI] = WEIGHTEDREGRESS(...,'ETYP',ETYP)
% Optionally, the error can be determined analytically or by bootstrapping
% (N=1000).
%
% See also REGSTATS, ROBUSTFIT, NORMPDF

% (c) 2011-04-26 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Check input
if size(X,2)>1
	X=X';
	if size(X,2)>1
		error('X should be a vector');
	end
end
mx = size(X,1);
[m,n] = size(Y); % n = number of indepedent variables
if m~=mx
	Y = Y';
	[m,n] = size(Y); %#ok<ASGLU> % n = number of indepedent variables
end
if size(Z,2)>1
	Z   = Z';
	if size(Z,2)>1
		error('Z should be a vector');
	end
end
if nargin<5 % Set a default
	XI = linspace(min(X),max(X),20);
end
if size(XI,2)>1
	XI=XI';
	if size(XI,2)>1
		error('X should be a vector');
	end
end
%% Optional, undocumented type of fit:
% simple regression (default) or a robust fitting procedure
fittype		= keyval('fittype',varargin,'regstats');
wfun		= keyval('wfun',varargin,'gaussian');
nbin		= keyval('nbin',varargin,0);
etyp		= keyval('etyp',varargin,'analytic');


switch wfun
	case 'gaussian'
		nreg	= n; % number of regression parameters. Bias is obsolete for gaussian weighting
		regindx = 2:(nreg+1);
	otherwise
		nreg	= n+1; % number of regression parameters.
		regindx = 1:nreg;
end

%% Weighted regress
% initialization of matrices
beta	= NaN(length(XI),nreg); % regression coefficients
se		= beta; % standard error
xi		= NaN(size(XI)); % wfuned average x-parameter
N = xi;
for ii   = 1:length(XI)
	switch wfun
		case 'gaussian'
			w1		= normpdf(X,XI(ii),sigma); % Gaussian wfuning for one column
			sw		= nansum(w1); % normalization factor
			w		= repmat(w1,1,n); % wfuning for number of independent measures
			xi(ii)	= dot(w1,X)./sw; % nansum(w1.*X)
			x		= w.*Y./sw;
			y		= w1.*Z./sw;
			sel		= X>=xi(ii)-2*sigma  & X<=xi(ii)+2*sigma;
		case 'boxcar'
			sel		= X>=XI(ii)-sigma  & X<=XI(ii)+sigma;
			xi(ii)	= nanmean(X(sel));
			x		= Y(sel,:);
			y		= Z(sel);
			N(ii) = sum(sel);
	end
	if nansum(sel)>nbin && nansum(sel)>n
		switch fittype
			case 'regstats'
				b	= regstats(y,x,'linear',{'beta';'tstat'});
				beta(ii,:)			= b.beta(regindx);
				
				switch etyp
					case 'boot'
						s = bootstrp(50, @regstats ,y,x); % bootstrap to obtain sd of sd of residuals
						s = [s.beta];
						s = std(s,[],2);
						se(ii,:) = s(regindx);
					otherwise
						se(ii,:)			= b.tstat.se(regindx);
				end
			case 'robustfit'
				[b,stats]	= robustfit(x,y);
				beta(ii,:)			= b(regindx);
				switch etyp
					case 'boot'
						s = std(bootstrp(50, @robustfit ,x,y)); % bootstrap to obtain sd of sd of residuals
						se(ii,:) = s(regindx);
					otherwise
						se(ii,:)			= stats.se(regindx);
				end
		end
	end
end

%% remove nans
% NaNs screw up interpolation
sel		= ~isnan(beta(:,1));
beta	= beta(sel,:);
se		= se(sel,:);
xi		= xi(sel);
N		= N(sel);
%% unique x-values
% interpolation dis not possible if there are multiple equivalent x-values
% solution: throw away equivalent x-values with corresponding y-values,
% keeping only the first value
[xi,indx] = unique(xi);
beta	= beta(indx,:);
se		= se(indx,:);
N		= N(indx);
whos N se
%% interpolate at XI

if numel(xi)>1 % to interpolate you should have at least 2 data points
	betai	= NaN(length(XI),size(beta,2));
	sei		= betai;
	for ii	= 1:size(beta,2)
		betai(:,ii) = interp1(xi,beta(:,ii),XI,'linear','extrap');
		sei(:,ii)	= interp1(xi,se(:,ii),XI,'linear','extrap');
	end
			Ni	= interp1(xi,N,XI,'linear','extrap');

else % Create some NaNs
	betai	= zeros(length(XI),size(beta,2));
	sei		= betai;
end
whos Ni sei


