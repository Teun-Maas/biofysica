function doublesound_MLR_prediction

%% Initialization
close all;
theta	= [0 0; 1 0; 0 0.5; 0.7 0.5];
ntheta	= size(theta,1);

sd		= 7;
x1		= -40:5:40;
x2		= -40:5:40;
[x1,x2] = meshgrid(x1,x2);

x1		= x1(:);
x2		= x2(:);
d		= x2-x1;

for ii	= 1:ntheta
	
	
	%% create
	Ravg	= binornd(1,theta(ii,1),size(x1));
	Rbin	= binornd(1,theta(ii,2),size(x1));
	y		= x1;
	sel		= logical(Ravg);
	y(sel)	= (x1(sel)+x2(sel))/2;
	sel		= ~sel & logical(Rbin);
	y(sel)	= x2(sel);
	y		= y+sd*randn(size(y));
	
	%% regress
	
	Y = y;
	X = [x1 x2]; % two single sounds (implicitly weighted average)
% 	X = [x1 (x1+x2)/2]; % single sound and unweighted average
% impossible, confound: X = [x1 x2 (x1+x2)/2]; % single sounds and
% unweighted average

	b = regstats(Y,X,'linear',{'beta','rsquare'});
	
	%% plot
	figure(2)
	subplot(ntheta,3,1+(ii-1)*3)
	plot(x1,y,'.');
	str = ['g_1 = ' num2str(b.beta(2),'%.2f')];
	text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
	
	subplot(ntheta,3,2+(ii-1)*3)
	plot(x2,y,'.');
	str = ['g_2 = ' num2str(b.beta(3),'%.2f')];
	text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
	
	
	
	subplot(ntheta,3,3+(ii-1)*3)
	plot((x1+x2)/2,y,'.');
	str = ['r^2 = ' num2str(b.rsquare,'%.2f')];
	text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
	
end

%%
titlestr = {'Single 1','Average','Bimodal','70% Avg + 30% Bi'};
for jj = 1:ntheta
	subplot(ntheta,3,3*(jj-1)+1)
	ylabel('response (deg)');
	
	subplot(ntheta,3,3*(jj-1)+2)
	title(titlestr{jj});

	for ii = 1:3
		subplot(ntheta,3,ii+(jj-1)*3)
		axis square
		xlim([-60 60]);
		ylim([-60 60]);
		box off
		% 	title(b.rsquare)
		set(gca,'XTick',-40:20:40,'YTick',-40:20:40);
		unityline('k:');
		horline(0,'k:');
	end
end

%% xlabel
subplot(ntheta,3,1+(ntheta-1)*3);
xlabel('stimulus 1 (deg)');

subplot(ntheta,3,2+(ntheta-1)*3);
xlabel('stimulus 2 (deg)');

subplot(ntheta,3,3+(ntheta-1)*3);
xlabel('average (deg)');

savegraph(mfilename,'eps');

% function writemodel
% 
% 
% model {
%     # Likelihood:
%     for( i in 1 : N ) {
%       y[i] ~ dnorm( mu[i] , tau ) 
%       mu[i] <- muOfClust[ clust[i] ]
%       clust[i] ~ dcat( pClust[1:Nclust] )
%     }
%     # Prior:
%     tau ~ dgamma( 0.01 , 0.01 )
%     for ( clustIdx in 1: Nclust ) {
%       muOfClust[clustIdx] ~ dnorm( 0 , 1.0E-10 )
%     }
%     pClust[1:Nclust] ~ ddirch( onesRepNclust )
% }
 
% The data specification:
% 
% # Generate random data from known parameter values:
% set.seed(47405)
% trueM1 = 100
% N1 = 200
% trueM2 = 145 # 145 for first example below; 130 for second example
% N2 = 200
% trueSD = 15
% effsz = abs( trueM2 - trueM1 ) / trueSD
% y1 = rnorm( N1 ) 
% y1 = (y1-mean(y1))/sd(y1) * trueSD + trueM1
% y2 = rnorm( N2 ) 
% y2 = (y2-mean(y2))/sd(y2) * trueSD + trueM2
% y = c( y1 , y2 ) 
% N = length(y)
% 
% # Must have at least one data point with fixed assignment 
% # to each cluster, otherwise some clusters will end up empty:
% Nclust = 2
% clust = rep(NA,N) 
% clust[which.min(y)]=1 # smallest value assigned to cluster 1
% clust[which.max(y)]=2 # highest value assigned to cluster 2 
% dataList = list(
%     y = y ,
%     N = N ,
%     Nclust = Nclust ,
%     clust = clust ,
%     onesRepNclust = rep(1,Nclust) 
% )