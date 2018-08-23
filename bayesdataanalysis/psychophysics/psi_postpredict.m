function psi_postpredict(x,y,s,samples,fun,centroidFlag,fig)
% PSI_POSTPREDICT(X,Y,S,SAMPLES,FUN,CENTROID,FIG)
%
% See also PSIFIT

%% Parameters
theta			= samples.theta;
omega			= samples.omega;
gamma			= samples.gamma;
lambda			= samples.lambda;

% guess		= samples.guess;

%% Determine centroid
if strcmp(centroidFlag,'mode')
	centroid = ['bf_' centroidFlag]; % for random samples, BF_MODE works better than MODE
else
	centroid = centroidFlag;
end
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	str = [parameterNames{parIdx} ' = samples.(parameterNames{parIdx});'];
	eval(str);
	str = ['c' parameterNames{parIdx} ' = ' centroid '(samples.(parameterNames{parIdx}));'];
	eval(str);
end

%%
chainLength = length(theta);

%% believable psychometric curve
cVec		= floor(linspace(1,chainLength,30));
xWid		= max(x)-min(x);
xComb		= linspace(min(x)-0.1*xWid,max(x)+0.1*xWid,201);
xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	selw		= xComb>=0;
	xComb		= xComb(selw);
end

%% some shenangins to plot a max of 25 subplots per figure
us = unique(s);
ns = length(us);
sb = ceil(sqrt(ns));
if sb>5
	sb = 5;
end
figcnt	= 0;
sbcnt	= 0;
cnt		= false;


for ii = 1:ns
	% for every subject/group
	sbcnt		= sbcnt+1; % new subplot count
	if cnt
		sbcnt	= 1;
		figcnt	= figcnt+1;
		cnt		= false;
	end
	if mod(ii,25)==0 % when 25 subplots have been reached
		cnt		= true; % new figure (see above)
	end
	sel			= s==us(ii);
	
	figure(fig+figcnt);
	subplot(sb,sb,sbcnt)
	hold on
	
	%% Posterior predictive
	% Data
	if size(y,2)==1 % y = bernouilli rate
		[r,ux]	= assemble(y(sel),x(sel),'fun',@sum);
		n		= assemble(ones(size(y(sel))),x(sel),'fun',@sum);
	elseif size(y,2)==2 % y = [rate n]
		 % % this should work
% 		r = y(sel,1);
% 		n = y(sel,2);
% 		ux = x(sel);
% % but perhaps the experimenter did not 'assemble' the data correctly
		[r,ux]	= assemble(y(sel,1),x(sel),'fun',@sum);
		n		= assemble(y(sel,2),x(sel),'fun',@sum);		
	end
	rprop			= r./n;
	[lb, ub]		= binomialci(r, n, 0.05);	
	%
% 	for cIdx	= cVec
% 		ypred	= psifun(ux,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
% 		whos n ypred
% 		keyboard
% 		R = binornd(n,ypred);
% 	end
	
	%% Posterior credible psychometric curves
	halfway = ((1-clambda)+cgamma)/2;
	cnt = 0;
	ypred = psifun(xComb,ctheta(ii),comega(ii),cgamma(ii),clambda(ii),0.1,'function',fun);
	M = NaN(numel(cVec),numel(xComb));
	for cIdx	= cVec
		cnt = cnt+1;
		ypred	= psifun(xComb,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
% 		xInt	= theta(cIdx,ii);

		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.9 .9 .9]);
% 		plot( [xInt xInt],[halfway(ii) -0.1], 'k:','Color',[.7 .7 .7]);
		
		M(cnt,:) = ypred;
	end

	mu = mean(M);
	sd = std(M);
	
	%% Max Posterior Curve
	ypred = psifun(xComb,ctheta(ii),comega(ii),cgamma(ii),clambda(ii),0.1,'function',fun);
	
	
	%% Graphics
	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
	errorpatch(xComb,mu,3*sd,'k');
	xInt		= ctheta(ii);
	plot([xInt xInt],[halfway(ii) -0.1], 'k-','Color','k');
	horline(cgamma(ii),'k:');
	horline(1-clambda(ii),'k:');
	str = {['\theta =' num2str(ctheta(ii),'%.1f') ', \omega= ' num2str(comega(ii),'%.2f')],...
		['\gamma = ' num2str(round(100*cgamma(ii))) '% , \lambda = ' num2str(round(100*clambda(ii))) '%']};
	text(mean(x),1.1,str,'HorizontalAlignment','center');
	
    
% 	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5); % data
	errorbar(ux,rprop,rprop-lb,ub-rprop,'ks','MarkerFaceColor','w','MarkerSize',10); % data
	
	hdi = hdimcmc(samples.theta);
	plot(hdi,[0.05 0.05],'k-','LineWidth',2);

	hdi = hdimcmc(samples.omega);
	plot([ctheta(ii)-hdi(2)/2 ctheta(ii)+hdi(2)/2],[-0.05 -0.05],'k-','LineWidth',2,'Color',[.7 .7 .7]);
	plot([ctheta(ii)-hdi(1)/2 ctheta(ii)+hdi(1)/2],[-0.05 -0.05],'k-','LineWidth',2);
	
% 	%% Density
% 	if size(y,2)==1
% 		if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
% 			supp = 'positive';
% 			xi		= linspace(0.1,max(x),100);
% 		else
% 			supp = 'unbounded';
% 			xi		= linspace(min(x),max(x),100);
% 		end
% 		% estimate rate from Bernouilli data via ksdensity
% 		sel		= y==0;
% 		f0		= ksdensity(x(sel),xi,'support',supp);
% 		
% 		sel		= y==1;
% 		f1		= ksdensity(x(sel),xi,'support',supp);
% 		
% 		f = f1-f0;
% 		f = f-min(f);
% 		f = f./max(f);
% 		f = (1-clambda(ii)-cgamma(ii))*f+cgamma(ii);
% 		
% 		plot(xi,f,'r-','LineWidth',2);
% 		
% 	end
	%% Labels
	title(us(ii))
	set(gca,'YTick',0:0.25:1);
% 	set(gca,'XTick',[ctheta-comega/2 ctheta ctheta+comega/2]);

	xlabel('X');
	ylabel('Probability');
	axis square;
	box off;
	ylim([-0.1 1.2]);
	xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
	nicegraph;
	grid off
end
