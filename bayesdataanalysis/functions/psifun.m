function p = psifun(x,theta,omega,gamma,lambda,alpha,varargin)
% P = PSIFUN(X,THETA,OMEGA,GAMMA,LAMBDA,ALPHA)
%
% Psychometric function of the form:
%
% p = (1-lambda)*( (1-gamma)*fun(x,theta,omega,alpha)+gamma ) + gamma*lambda;
%
% See also PSIFIT, LOGISTICFUN, PROBITFUN, GUMBELFUN, REVGUMBELFUN,
% WEIBULLFUN, REVWEIBULLFUN

%% Initialization
hfun = keyval('function',varargin,@logisticfun);

if ~isa(hfun,'function_handle')
	error('Please pass a function handle')
end
strfun	= func2str(hfun);
sel		= strcmp(strfun(1:end-3),{'logistic','probit','gumbel','revgumbel','weibull','revweibull'});


if any(sel)
	p = (1-lambda)*( (1-gamma)*hfun(x,theta,omega,alpha)+gamma ) + gamma*lambda;
else 
	error([upper(strfun) ' does not exist: choose another psi function']);
end
