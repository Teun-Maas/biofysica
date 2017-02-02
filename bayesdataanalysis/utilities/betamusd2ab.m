function betaAB = betamusd2ab(mu,sd)
% betaAB = BETAMUSD2AB(MU,SD)
%
% Obtain the beta shape parameters betaAB.A and betaAB.B from the mean MU and standard
% deviation SD

% Original in R:	Kruschke, J. K. (2014). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij

if mu<=0 || mu>1
	msgid		= 'bayesian:betaadfrommeansd:meanOutOfRange';
	errstr		= 'Must have 0 < mean < 1';
	error(msgid,errstr);
end
if sd<0 
	msgid		= 'bayesian:betaadfrommeansd:stdOutOfRange';
	errstr		= 'Standar deviation must be >0';
	error(msgid,errstr);
end

kappa = mu.*(1-mu)./sd.^2 - 1;
if kappa<=0
	msgid		= 'bayesian:betaadfrommeansd:invalidmeanandstd';
	errstr		= 'invalid combination of mean and sd';
	error(msgid,errstr);
end
betaAB.a	= mu.*kappa;
betaAB.b	= (1-mu).*kappa;
