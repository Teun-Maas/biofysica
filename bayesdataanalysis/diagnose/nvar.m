function [nVar,varNames] = nvar(mcmc)
% N = NVAR(X)
%
% gives number of variables of an MCMC struct
%
% [N,VARNAME] = NVAR(X)
%
% also outputs variable names VARNAME if X is a structure (inputname otherwise)
%
% See also NCHAIN, NITER


if isstruct(mcmc)
	varNames	= fieldnames(mcmc);
	nVar		= numel(varNames);
else
	nVar		= 1;
	varNames	= [];
end

