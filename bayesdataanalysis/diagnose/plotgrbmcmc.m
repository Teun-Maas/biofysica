function plotgrbmcmc(mcmcStruct,varargin)
% PLOTGRBMCMC(MCMC)
%
% Gelman-Rubin-Brooks plot
% shows the evolution of Gelman and Rubin's shrink factor as the number of
% iterations increases. 
%
%The ?potential scale reduction factor? is calculated for each variable in
%x, together with upper and lower confidence limits. Approximate
%convergence is diagnosed when the upper limit is close to 1. For
%multivariate chains, a multivariate value is calculated that bounds above
%the potential scale reduction factor for any linear combination of the
%(possibly transformed) variables.
% The confidence limits are based on the assumption that the stationary
% distribution of the variable under examination is normal. Hence the
% ?transform? parameter may be used to improve the normal approximation.    
%
%
% References:
% Gelman, A and Rubin, DB (1992) Inference from iterative simulation using
% multiple sequences, Statistical Science, 7, 457-511. 
%
% Brooks, S P. and Gelman, A. (1998) General Methods for Monitoring
% Convergence of Iterative Simulations. Journal of Computational and
% Graphical Statistics, 7, 434-455.  
%
%
% Above text comes from gelman.plot and gelman.diag in R
%
% See also DIAGMCMC, SHRINKFACTOR

%% Initialization
parNames	= fieldnames(mcmcStruct);
parName		=  keyval('parName',varargin,parNames{1});


[rhat,rhat95] = shrinkfactor(mcmcStruct,'parName',parName,'binwidth',1);
t = 1:length(rhat);
%% Graphics
plot(t,rhat,'k-');
hold on
plot(t,rhat95,'k:');

title(parName);

box off;
xlim([-100 max(t)+100]);
horline(1,'k-');
xlabel('Last iteration in chain');
ylabel('Shrink factor')
legend({'median','97.5%'});