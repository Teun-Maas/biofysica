function x = invpsifun(p,theta,omega,gamma,lambda,alpha,varargin)
% X = INVPSIFUN(P,THETA,OMEGA,GAMMA,LAMBDA,ALPHA)
%
% Inverse of the psychometric function of the form:
%
% p = (1-lambda)*( (1-gamma)*fun(x,theta,omega,alpha)+gamma ) + gamma*lambda;
%
% Currently only works for the logistic function
% See also PSIFIT, LOGISTICFUN

%% Initialization
hfun = keyval('function',varargin,@logisticfun);

if ~isa(hfun,'function_handle');
	error('Please pass a function handle')
end
strfun	= func2str(hfun);
% sel		= strcmp(strfun(1:end-3),{'logistic','probit','gumbel','revgumbel','weibull','revweibull'});

sel		= strcmp(strfun(1:end-3),{'logistic'});

if any(sel)
	% 	hfun(x,...) = 1./(1+exp(-z/omega*(x-theta)));
	z = 2*log(1/alpha-1);
	% 	p = (1-lambda)*( (1-gamma)*hfun(x,theta,omega,alpha)+gamma ) + gamma*lambda;
	% 	(1-lambda)*( (1-gamma)*hfun(x,theta,omega,alpha)+gamma ) = p - gamma*lambda;
	% 	(1-gamma)*hfun(x,theta,omega,alpha)+gamma = (p - gamma*lambda)/(1-lambda);
	% 	(1-gamma)*hfun(x,theta,omega,alpha) = (p - gamma*lambda)/(1-lambda) - gamma;
	% 	hfun(x,theta,omega,alpha) = ((p - gamma*lambda)/(1-lambda) - gamma)/(1-gamma);
	% 	1./(1+exp(-z/omega*(x-theta))) = ((p - gamma*lambda)/(1-lambda) - gamma)/(1-gamma);
	% 	1+exp(-z/omega*(x-theta)) = (1-gamma)/((p - gamma*lambda)/(1-lambda) - gamma);
	% 	exp(-z/omega*(x-theta)) = (1-gamma)/((p - gamma*lambda)/(1-lambda) - gamma)-1;
	% 	-z/omega*(x-theta) = log((1-gamma)/((p - gamma*lambda)/(1-lambda) - gamma)-1);
	% 	x-theta = log((1-gamma)/((p - gamma*lambda)/(1-lambda) - gamma)-1)*omega/-z;
		x = log((1-gamma)./((p - gamma*lambda)./(1-lambda) - gamma)-1)*omega./-z+theta;
% 	inlog	= ((1-lambda)*(1-gamma)-(p-gamma))./(p-gamma);
% 	x		= -omega./z*log(inlog)+theta;
% 	x= abs(x);
else
	error([upper(strfun) ' does not exist: choose another invpsi function']);
end
