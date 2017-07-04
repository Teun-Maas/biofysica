function [mui,sei,XI] = weightedmean(X,Y,sigma,XI,varargin)
% [mu,SE,XI] = PA_WEIGHTEDMEAN(X,Y,SIGMA,XI)
%
% Perform linear regression  Y vs Z weighted by factor X.
% Weighing is performed either through a Gaussian or boxcar profile with a
% width of SIGMA, centered  around XI.
%
% [mu,SE,XI] = PA_WEIGHTEDREGRESS(...,'NBIN',NBIN)
% Optionally, you can provide the minimum number NBIN of data points required
% within plusmin 1 SIGMA of XI (default = 0). This is useful when few
% data-points ara available at XI, and regression might be biased towards
% outliers.
%
% [mu,SE,XI] = PA_WEIGHTEDREGRESS(...,'WFUN',WFUN)
% Optionally, you can provide weight function: 'gaussian' or 'boxcar'
% (default WFUN = 'gaussian').
%
%
% See also REGSTATS, ROBUSTFIT, NORMPDF

% (c) 2011-04-26 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

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

if nargin<4 % Set a default
	XI = linspace(min(X),max(X),20);
end
if size(XI,2)>1
	XI=XI';
	if size(XI,2)>1
		error('X should be a vector');
	end
end
sel = ~isnan(Y);
Y = Y(sel);
X = X(sel);

%% Optional, undocumented type of fit:
% simple regression (default) or a robust fitting procedure
wfun         = keyval('wfun',varargin);
if isempty(wfun)
	wfun = 'gaussian';
end
nbin        = keyval('nbin',varargin);
if isempty(nbin)
	nbin = 0;
end
nboot         = keyval('nboot',varargin);

%% Weighted regress
switch wfun
	case 'gaussian'
		[mu,se,xi] = getmean(X,Y,XI,sigma,nbin,n);
		% 		se		= getsd(Y,X,xi,sigma);
		%         CI		=  NaN(length(XI),2);
		%         for ii   = 1:length(XI)
		%             ci			= bootstrp(nboot, @getmean ,X',Y',XI(ii),sigma,nbin,n);
		% %             CI(ii,:) = bootci(nboot,@getmean ,X,Y,XI(ii),sigma,nbin,n)
		%             CI(ii,:) = prctile(ci,[45 55]);
		% %             keyboard
		%         end
		%         se = CI;
	case 'boxcar'
		% initialization of matrices
		mu	= NaN(length(XI),1); % regression coefficients
		se		= mu; % standard error
		xi		= NaN(size(XI)); % wfuned average x-parameter
		for ii   = 1:length(XI)
			sel		= X>=XI(ii)-sigma  & X<=XI(ii)+sigma;
			xi(ii)	= nanmean(X(sel));
			y		= Y(sel);
			mu(ii) = nanmean(y);
			
			se(ii) = nanstd(y);
		end
end

%% remove nans
% NaNs screw up interpolation
% sel		= ~isnan(mu(:,1));
% mu      = mu(sel,:);
% se		= se(sel,:);
% xi		= xi(sel);

%% unique x-values
% interpolation dis not possible if there are multiple equivalent x-values
% solution: throw away equivalent x-values with corresponding y-values,
% keeping only the first value
[~,indx] = unique(xi);
mui      = mu(indx,:);
sei		= se(indx,:);
XI = XI(indx);
% %% interpolate at XI
% if numel(xi)>1 % to interpolate you should have at least 2 data points
% 	mui	= NaN(length(XI),size(mu,2));
% 	sei	= NaN(length(XI),2);
% 	mui = interp1(xi,mu,XI,'linear');
% 	if ~strcmpi(wfun,'boxcar');
% 		sei(:,1)	= interp1(xi,se(:,1),XI,'linear');
% 		sei(:,2)	= interp1(xi,se(:,2),XI,'linear');
% 	else
% 		sei = se;
% 	end
% else % Create some NaNs
% 	xi = XI;
% 	mui	= zeros(length(XI),size(mu,2));
% 	sei		= mui;
% 	
% end

function [mu,se,xi] = getmean(X,Y,XI,sigma,nbin,n)
% initialization of matrices
mu      = NaN(length(XI),1); % regression coefficients
se		= mu; % standard error
xi		= NaN(size(XI)); % wfuned average x-parameter
for ii   = 1:length(XI)
	w1		= normpdf(X,XI(ii),sigma); % Gaussian wfuning for one column
	sw		= nansum(w1); % normalization factor
	xi(ii)	= nansum(w1.*X)./sw;
	xi(ii) = XI(ii);
	y		= w1.*Y./sw;
	sel		= X>=XI(ii)-sigma  & X<=xi(ii)+sigma;
	if sum(sel)>nbin && sum(sel)>n
		mu(ii) = nansum(y);
	end
	
	se(ii) = getsd(Y,X,xi(ii),sigma);
end


function S = getsd(Y,X,XI,sigma)

w	= normpdf(X,XI,sigma);
sw	= sum(w);
V	= 0;
for k	= 1:length(Y)
	V	= V+w(k)*(Y(k))^2;
end
n = sqrt(sum(numel(X).*w))/1.96;
S		= sqrt(V/sw)./n;