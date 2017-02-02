close all
x		= 0:50;
T		= 6;
lambdaA	= 6.6697;
lambdaV	= 3.9;

% T		= 0;
% lambdaA	= 0.8;
% lambdaV	= 0.1;


lambdaAV	= lambdaA+lambdaV;

%% Auditory



Y		= poisspdf(x,lambdaA);
C		= poisscdf(x,lambdaA);

sel		= x>T;
PA		= sum(Y(sel));

subplot(231)
stem(x,Y,'o-');
hold on
axis square;
box off
set(gca,'TickDir','out');
title(PA)
xlabel('# cues detected');
ylabel('Prob(#)');

subplot(212)
plot(x,C,'-');
hold on
axis square;
box off
set(gca,'TickDir','out');
title(PA)


%% Visual

Y		= poisspdf(x,lambdaV);
C		= poisscdf(x,lambdaV);
sel		= x>T;
PV		= sum(Y(sel));

subplot(232)
stem(x,Y,'ro-');
hold on
axis square;
box off
set(gca,'TickDir','out');
title(PV)
xlabel('# cues detected');
ylabel('Prob(#)');

subplot(212)
plot(x,C,'r-');

axis square;
box off
set(gca,'TickDir','out');
title(PV)

%% AudioVisual

Y		= poisspdf(x,lambdaAV);
C		= poisscdf(x,lambdaAV);
sel		= x>T;
PAV		= sum(Y(sel));

subplot(233)
stem(x,Y,'go-');
hold on
axis square;
box off
set(gca,'TickDir','out');
title(PAV)
xlabel('# cues detected');
ylabel('Prob(#)');

subplot(212)
plot(x,C,'g-');

axis square;
box off
set(gca,'TickDir','out');
title(PAV)


%%
Pprobcomb =  PA+PV-PA*PV
