close all
clearvars

% xx xx | HAx HAx | HAx HACI | xx xxCI
gain	= [.6 .4 .6 .7];
bias	= [0 -50 -30 +50];
sd		= [12 12 14 12];

nsubjects	= 9;
ntrials		= 45;
nlevels		= 1;
x			= -70:10:70;
x			= repmat(x,1,nlevels);
nexp		= 1000;
ncond		= 4;
P = NaN(nexp,1);
for ii = 1:nexp
	MAE = NaN(ncond,nsubjects);

	for jj = 1:ncond
		g = gain(jj)+0.1*(2*rand(nsubjects,1)-1);
		b = bias(jj)+10*(2*rand(nsubjects,1)-1);
		s = sd(jj)+3*(2*rand(nsubjects,1)-1);
		
		for kk = 1:nsubjects
		y = g(kk)*x+b(kk)+s(kk)*randn(size(x));
		MAE(jj,kk) = mean(abs(y-x));
% 		b = regstats(y,x,'linear','beta');
		end
	end
	
	[H,P(ii)] = ttest(MAE(1,:),MAE(2,:));
end

%%

figure(1)
clf
plotpost(P,'showCurve',true)

title(round(100*sum(P<0.05)/nexp))