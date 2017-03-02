function pa_how2race
% PA_HOW2RACE
% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

close all
clear all
clc

n = 1000;
rt1 = sort(200+50*randn(n,1));
rt2 = sort(200+70*randn(n,1));
rt3 = sort(150+40*randn(n,1));

p = 0.1:0.1:0.9;
q1 = quantile(rt1,p);
q2 = quantile(rt2,p);
q3 = quantile(rt3,p);

P = (1:n)/n;
plot(rt1,P,'-','Color',[.7 .7 .7]);
hold on
plot(rt2,P,'-','Color',[.7 .7 .7]);
plot(rt3,P,'k-');

plot(q1,p,'ko','MarkerFaceColor',[.7 .7 .7],'LineWidth',2);
plot(q2,p,'ko','MarkerFaceColor',[.7 .7 .7],'LineWidth',2);
plot(q3,p,'ko','MarkerFaceColor','r','LineWidth',2);

t = 0:1:500;
nt = numel(t);
c1 = NaN(nt,1);
c2 = c1;
c3 = c1;
for ii = 1:nt
	c1(ii) = sum(rt1<t(ii));
	c2(ii) = sum(rt2<t(ii));
	c3(ii) = sum(rt3<t(ii));
end
c1 = c1/max(c1);
c2 = c2/max(c2);
c3 = c3/max(c3);

plot(t,c1,'.','Color',[.7 .7 .7]);
plot(t,c2,'.','Color',[.7 .7 .7]);
plot(t,c3,'.','Color','k');

ks1 = ksdensity(rt1,t,'function','cdf');
ks2 = ksdensity(rt2,t,'function','cdf');
ks3 = ksdensity(rt3,t,'function','cdf');

plot(t,ks1,'-','LineWidth',2,'Color',[.7 .7 .7]);
plot(t,ks2,'-','LineWidth',2,'Color',[.7 .7 .7]);
plot(t,ks3,'k-','LineWidth',2);

ulrich = ks1+ks2;
sel = ulrich>1;
ulrich(sel) = 1;
plot(t,ulrich,'k-','LineWidth',2);

gielen = ks1+ks2-ks1.*ks2;
plot(t,gielen,'k.-','LineWidth',2);

colonius = max(ks1,ks2);
plot(t,colonius,'k-','LineWidth',2);
X = t;
Y = gielen;
E = [colonius; ulrich];
whos X Y E
pa_errorpatch(X,Y,E,'r');

xlim([0 500]);
axis square;
box off
xlabel('Reaction time (ms)');
ylabel('Cumulative probability P(\tau<t)');