function double_regress

%% Initialization
close all

%% The data specification:
delay	= [0 5 10 20 40 80 120 160 320]';
ndelays = numel(delay);
subjects = 1:9;
nsubj = numel(subjects);
G = NaN(nsubj,ndelays);

SD = G;
for jj = 1:nsubj
	for ii = 1:ndelays
		dT		= delay(ii);
		[r,x,y1,x1] = getdata(subjects(jj),dT);
		
		
		% single linear regression
		% r1 = g*x1
		b1			= regstats(y1,x1,'linear','beta');
		% 		b1.beta		= [0 1];
		% 				b1.beta(2)		=  1;
		
		% 'multiple' linear regression
		% r = a*p1+(1-a)*p2+b
		% r = a*p1+p2-a*p2+b
		% r-p2 = a*(p1-p2)+b
		% r-g*x2 = a*g*(x1-x2)+b
		% 'y = a*x+b'
		
		p1 = b1.beta(1)+b1.beta(2)*x(:,1);
		p2 = b1.beta(1)+b1.beta(2)*x(:,2);
		
		
		
		y = r-p2;
		x = p1-p2;
		b = regstats(y,x,'linear',{'beta','tstat','r'});
		G(jj,ii) = b.beta(2);
		Gse(jj,ii) = b.tstat.se(2);
		SD(jj,ii) = std(b.r);
		
		
	end
	
	G1(jj) = b1.beta(2);
	
	figure(1)
	subplot(3,3,jj)
	mu		= G(jj,:);
	sd		= Gse(jj,:);
	x		= 1:ndelays;
	hold on
	errorpatch(x,mu,sd)
	plot(x,mu,'ko','MarkerFaceColor','w');
	axis square;
	xlim([0 10]);
	ylim([-0.1 1.1]);
	horline([0 0.5 1]);
	horline(G1(jj),'r-');
	box off;
	set(gca,'TickDir','out','XTick',x,'XTickLabel',delay);
	xlabel('Delay (ms)');
	ylabel('Weight');
	
	% 	savegraph([mfilename num2str(jj)],'eps');
end
savegraph([mfilename '1']);

%%


figure(100)
clf
subplot(131)
mu		= mean(G);
sd		= std(G);
x		= 1:ndelays;
hold on
errorpatch(x,mu,sd)
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
xlim([0 10]);
ylim([-0.1 1.1]);
horline([0 0.5 1]);
box off;
set(gca,'TickDir','out','XTick',x,'XTickLabel',delay);
xlabel('Delay (ms)');
ylabel('Weight');


subplot(132)
mu		= mean(SD);
sd		= std(SD);
x		= 1:ndelays;
hold on
errorpatch(x,mu,sd)
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
xlim([0 10]);
ylim([-5 30]);
horline(0);
box off;
set(gca,'TickDir','out','XTick',x,'XTickLabel',delay);
xlabel('Delay (ms)');
ylabel('\sigma (deg)');

subplot(133)
plot(mean(SD),mean(G),'ko');
axis square
box off
set(gca,'TickDir','out');
xlabel('\sigma (deg)');
ylabel('Weight');
ylim([-0.1 1.1]);
xlim([0 20])
%%
savegraph(mfilename);
%% Stats
% kruskal-Wallis (non-parametric 'ANOVA')

[p,tbl,stats] = kruskalwallis(G);
%multcompare(stats)

[p,tbl,stats] = kruskalwallis(SD);
%multcompare(stats)


function [y,x,y1,x1] = getdata(subject,dT)
load('/Users/marcw/Dropbox/Work/manuscripts/Rachel PhD/Papers/Paper 2 A9 Double Sounds/Data/A9DataMatrix.mat');


% 1) Trialnumber
% 2) 1snd location (el)
% 3) 2snd location (el)
% 4) el response
% 5) az resp (az location is always zero so not listed in matrix)
% 6) 1.snd
% 7) 2.snd
% 8) subjID (1-9)
% 9) delay (1:17 (-320:320)),
% 10) 1=single sound trials 2= bzz vs gwn trials 3= bzz&gwn (same location)
delay	= -[-320 -160 -120 -80 -40 -20 -10 -5 0 5 10 20 40 80 120 160 320]';


sel			= ismember(delay(Matrix(:,9)),[dT -dT]) & Matrix(:,10)==2;
y			= Matrix(sel,4);
dt			= delay(Matrix(sel,9));
s			= Matrix(sel,8);
RT			= Matrix(sel,11);
x1			= Matrix(sel,2);
x2			= Matrix(sel,3);
x			= [x1 x2];
sel			= dt==-dT;
x(sel,:)	= [x2(sel) x1(sel)];


sel			= dt>=0;
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);
sum(sel)

sel			= ~isnan(RT);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

sel			= ismember(s,subject);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

d			= abs(x(:,2)-x(:,1));
sel			= d>20;
x			= x(sel,:);
y			= y(sel);
RT			= RT(sel);

hdi			= hdimcmc(RT);
sel			= RT>hdi(1) & RT<hdi(2);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

N			= numel(y);

%% Single
sel		= ismember(Matrix(:,10),[1 3]);
y1		= Matrix(sel,4);
x1		= Matrix(sel,2);
N1		= numel(y);
s1		= Matrix(sel,8);