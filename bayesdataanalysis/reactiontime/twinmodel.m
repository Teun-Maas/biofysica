function tmp
% Visual stimulus is target, auditory is distracter
% t = SoA, o = oidth of integration oindoo
% t < t+o <0
clear all
close all

%% Focused attention oithout oarning
T = -500:500;
% T = -T;
d = 142;
m = 193;
o = 342;
lv = 1/119;
la = 1/45;
ga = 300;
k = 20;
Par(1) = lv;
Par(2)= la;
Par(3) = m;
Par(4) = o;
Par(5) = d;
Par(6) = ga;
Par(7) = k;


%%

figure(1)
ax(1) = subplot(131);
PI = twinpi(T,lv,la,o);
RT = 1./lv+m-d*PI;
plot(T,RT)

axis square

ax(2) = subplot(132);
PW = twinpw(T,lv,la,ga);
RT = 1./lv+m-k*PW;
plot(T,RT)
axis square

ax(3) = subplot(133);
PW = twinpw(T,lv,la,ga);
RT = 1./lv+m-d*PI-k*PW;
plot(T,RT)
axis square

linkaxes(ax);

for ii = 1:3
	subplot(1,3,ii)
	horline(1/lv+m);
end
axis auto

RT = twinfun(T,Par);

% RT = 1./lv+m-d*PI-k*PW;
plot(T,RT,'r-')
axis square

%% twinfit
T = [-600 -300:100:200 600];
% T = -500:60:500;

RT = twinfun(T,Par);
% RT = 1./lv+m-d*PI-k*PW;


% RT1 = 1./lv+m;
% RT2 = 1./la+m;
% 
% RT1 = 1./lv+m+10*randn;
% RT2 = 1./la+m+10*randn;
RT = RT+20*randn(size(RT));

hold on
plot(T,RT,'ro')

drawnow
%%
% twinfit(T,RT,RT1,RT2)
P = fittwin(T,RT);

RT = twinfun(T,Par);
plot(T,RT,'ko')

% keyboard
function Par = fittwin(X,Y,x0)
%FITSIN   Fit a sine through data.
%   [PAR]=FITSIN(X,Y) returns parameters of the sine a*sin(b*X+c)+d,
%   in the following order: [a b c d].
%
%   See also SINFUN, SINERR.
%
%  drs M.
if nargin<3
x0(1) = 1/119;
x0(2) = 1/45;
x0(3) = 193;
x0(4) = 342;
x0(5) = 142;
x0(6) = 300;
x0(7) = 20;
x0 = [1/150 1/150 0 0 0 0 0];
LB = [1/150 1/150 0 0 0 0 0];
UB = [1/5 1/5 3000 600 Inf Inf Inf];
% 	parstart = [lv		la		m	o	d	ga	k;
% 				1/150	1/150	0	0	0	0	0;
% 				1/5		1/5		Inf 600 Inf Inf Inf];
end
try
Par=fminsearchbnd(@twinerr,x0,LB,UB,[],X,Y);
catch
	Par=fminsearch(@twinerr,x0,[],X,Y);
	disp('Uh-oh');

end
% function y = sinfun(X,Par)
% % SINFUN   Function of the sine a*sin(b*X+c)+d.
% a = Par(1);
% b = Par(2);
% c = Par(3);
% d = Par(4);
% y = a*sin(b*X+c)+d;

function y = twinfun(X,Par)
% RT = TWINFUN(T,LAMBDA1,LAMBDA2,OMEGA,DELTA,GAMMA,KAPPA)
%
% See also TWINPI, TWINPW

lambda1 = Par(1);
lambda2 = Par(2);
mu		= Par(3);
omega	= Par(4);
delta	= Par(5);
gamma	= Par(6);
kappa	= Par(7);

% Probability of integration
PI = twinpi(X,lambda1,lambda2,omega);
% Probability of warning
PW = twinpw(X,lambda1,lambda2,gamma);
% Observed mean reaction time
y = 1./lambda1+mu-delta*PI-kappa*PW;

function err =  twinerr(Par,X,Y)
%SINERR   Determines error between experimental data and calculated sine.
%   [ERR]=SINERR(PAR,X,Y) returns the error between the calculated parameters
%   PAR, given by FITSIN and the parameters given by experimental data X and Y.
err = norm(Y-twinfun(X,Par));
