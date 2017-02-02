% Visual stimulus is target, auditory is distracter
% t = SoA, o = oidth of integration oindoo
% t < t+o <0
close all
clear all
%% Focused attention oithout oarning
T = -1000:10:1000;
% T = -T;
d = 142;
m = 193;
o = 342;
lv = 1/119;
la = 1/45;
ga = 300;
k = 20;


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

% for ii = 1:3
% 	subplot(1,3,ii)
% 	horline(1/lv+m);
% end
axis auto

RT = twinfun(T,lv,la,m,o,d,ga,k);
% RT = 1./lv+m-d*PI-k*PW;
plot(T,RT,'r-')
axis square

drawnow

[yfit,stats]	=	fittwin(T,RT);

% myfig;
figure;
plot(T,RT,'r.')
hold on
plot(T,yfit,'k-')
axis('square')
