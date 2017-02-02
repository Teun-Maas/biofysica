function [yfit,varargout] = fittwin(xin,yin,parstart)
%
%	[yfit,stats] = fitexp(xin,yin,[type],[parstart])
%	
%	See also fminsearch
%
%	INPUT:
%	xin		 =	x values
%	yin		 =	y values
%	type	 =	0: yfit = par(1) .* exp(-X .* par(3)) + par(2);
%				1: yfit = par(1) .* (X .^ par(3)) + par(2);
% 				default: 0
%	parstart =	[x0 x0 x0;
%				lb lb lb;
%				ub ub ub];
%				default:	[   0     0     0;
% 							-1000 -1000 -1000;
% 							 1000  1000  1000]
%
%	OUTPUT:
%	yfit	 =	depends on type (see above);
%	stats	 =	Optional. Returns a structure containing 
%				  R^2, p-value, SSR, SSE, N and coefficients in 'Coef'
%				  
%
%	Based on a function written by JohnvO
%	Modified for linear fits on 04.Nov 2015 by PeterBr

global X Y

X		=	xin;
Y		=	yin;

if nargin < 3

% For the estimation routine, ?V and ?A were restricted to a range consistent
% with neurophysiological estimates for peripheral processing times (Stein and Meredith 1993; Groh and Sparks, 1996)
% 5 ? 1/? ? 150.
% The width of the time window of integration,was restricted to a positive real number with an upper bound of 600 ms.
% Mean time for the second stage, ?, was restricted to be positive. 

% lambda1,lambda2,mu,omega,delta,gamma,kappa
d	=	100;%142;
m	=	100;%193;
o	=	100;%342;
lv	=	1/150;%1/119;
la	=	1/150;%1/45;
ga	=	100;%300;
k	=	10;%20;

% d	=	0;%142;
% m	=	0;%193;
% o	=	0;%342;
% lv	=	0;
% la	=	0;
% ga	=	0;%300;
% k	=	0;%20;
	parstart = [lv		la		m	o	d	ga	k;
				1/150	1/150	0	0	0	0	0;
				1/5		1/5		Inf 600 Inf Inf Inf];
end

% fitfcn = @(p,xdata)p(1) .* exp(-xdata .* p(3))+p(2);
% problem				=	createOptimProblem('lsqcurvefit','objective',@fitfcn, ...
% 											'x0',parstart(1,:),'lb',parstart(2,:),'ub',parstart(3,:), ...
% 											'xdata',xin,'ydata',yin);
										
%-- Global search for optimal parameters --%
problem				=	createOptimProblem('fmincon','objective',@getmse, ...
											'x0',parstart(1,:),'lb',parstart(2,:),'ub',parstart(3,:));

gs					=	GlobalSearch('Display','off');
[par,~,exitflag,output]	=	run(gs,problem);

%-- Calculate fit --%
yfit	=	twinfun(X,par(1),par(2),par(3),par(4),par(5),par(6),par(7));

%-- Calculation: Correlation Coefficient --%
%-- see http://mathworld.wolfram.com/CorrelationCoefficient.html --%
SSR		=	sum ((yfit - mean(Y)).^2);		%-- Sum of squared residuals --%
SSE		=	sum ((yfit - Y).^2);			%-- Sum of squared errors --%
R2		=	(SSR/(SSE+SSR))^(1/2);			%-- Correlation Coefficient R^2 --%

[~,P]	=	corrcoef(Y,yfit);

stats = struct('R2',R2,'p',P(2,1),'SSR',SSR,'SSE',SSE,'N',length(Y),'Coef',par);
varargout(1) = {stats};

if( exitflag ~= 1 )
	disp(output.message)
end




return

%-- Calculate mse for, in this case, an exponential --%
function mse = getmse(p)

global X Y

lambda1	=	p(1);
lambda2	=	p(2);
mu		=	p(3);
omega	=	p(4);
delta	=	p(5);
gamma	=	p(6);
kappa	=	p(7);

y		=	twinfun(X,lambda1,lambda2,mu,omega,delta,gamma,kappa);

mse		=	sum( (y-Y).^2); % squared error
   
return

function RT = twinfun(T,lambda1,lambda2,mu,omega,delta,gamma,kappa)
% RT = TWINFUN(T,LAMBDA1,LAMBDA2,OMEGA,DELTA,GAMMA,KAPPA)
%
% See also TWINPI, TWINPW

% Probability of integration
PI = twinpi(T,lambda1,lambda2,omega);
% Probability of warning
PW = twinpw(T,lambda1,lambda2,gamma);
% Observed mean reaction time
RT = 1./lambda1+mu-delta*PI-kappa*PW;

function PI = twinpi(T,l1,l2,o,varargin)
% PI = TWINPI(T,L1,L2,O)
%
% Probability of integration PI of the time-window-of-integration (TWIN)
% model, given stimulus onset asynchrony T, rate of the target L1, and rate
% of the distractor (in focused attention tasks) or target (in redundant
% target tasks) L2, and width of the time window O.
%
% For model description see e.g. ???
%
% See also TWINPW

task = keyval('task',varargin,'FA'); % task = focused attention (FA) or redundant target (RT)
n	= numel(T);
PI	= NaN(n,1);

switch task
	case 'FA'
		for tIdx = 1:n
			t = T(tIdx);
			if t<(t+o) && (t+o)<0
				PI(tIdx) = l1./(l1+l2).*exp(l2.*t).*(-1+exp(l2.*o));
			elseif t<=0 && 0<=(t+o)
				PI(tIdx) = 1./(l1+l2).*(l2*(1-exp(-l1.*(o+t)))+l1*(1-exp(l2.*t)));
			elseif 0<t && t<(t+o)
				PI(tIdx) = l2./(l1+l2).*(exp(-l1*t)-exp(-l1*(o+t)));
			end
		end
	case 'RT'
		disp('Not yet implemented');
end

function PW = twinpw(T,l1,l2,ga,varargin)
% PW = TWINPW(T,L1,L2,GA)
%
% Probability of warning PW of the time-window-of-integration (TWIN)
% model, given stimulus onset asynchrony T, rate of the target L1, and rate
% of the distractor (in focused attention tasks) or target (in redundant
% target tasks) L2, and GA.
%
% For model description see e.g. ???
%
% See also TWINPI

task = keyval('task',varargin,'FA'); % task = focused attention (FA) or redundant target (RT)
n	= numel(T);
PW	= NaN(n,1);


switch task
	case 'FA'
		for tIdx = 1:n
			t = T(tIdx);
			if (t+ga)<0
				PW(tIdx) = 1-(l1./(l1+l2).*exp(l2*(t+ga)));
			elseif (t+ga)>=0
				PW(tIdx) = (l2./(l1+l2).*exp(-l1*(t+ga)));
			end
		end
	case 'RT'
		disp('Not yet implemented');
end
