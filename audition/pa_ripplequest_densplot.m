%% Initialization
close all
clear all
clc

%% Load the data
pa_datadir;
cd('Test');
fname = 'ripplequest_EGJ2';
fname = 'ripplequest-MW-2012-06-21-2';
fname = 'ripple';

load(fname);

%% Determine threshold based on responses
n		= numel(Q);
vel		= NaN(n,1);
dens	= NaN(n,1);
mod		= NaN(n,1);
modsd = mod;
R = [];
L = [];
I = [];
for ii = 1:n
	q = Q(ii).q;

	% Ask Quest for the final estimate of threshold.
	t				= QuestMean(q);		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
	sd				= QuestSd(q);
	lat				= 1000*q.lat;
	try
		ntrials			= q.trialCount;
	catch
		ntrials = length(q.intensity);
	end
	response	= q.response(1:ntrials);
	intensity	= 10.^(q.intensity(1:ntrials));
	I			= [I; intensity]; %#ok<*AGROW> % modulation depth
	R			= [R; response]; % response
	L			= [L; lat]; % reaction time
	
	
	mod(ii)		= 10^t; % convert from logaritmic to linear
	modsd(ii)	= 10^sd;
	vel(ii)		= Q(ii).vel;
	dens(ii)	= Q(ii).dens;
	fprintf('Final threshold estimate (mean±sd) is %.2f ± %.2f\n',mod(ii),modsd(ii));
	
end
[dens,indx] = sort(dens);
mod			= -mod(indx)/100;
modsd		= modsd(indx)/100;
vel			= vel(indx);
I			= I(indx,:);
L			= L(indx,:);
R			= R(indx,:);

%% Plot threshold as a function of ripple density
figure(1)
errorbar(dens,mod,modsd,'ko-','MarkerFaceColor','w','LineWidth',2);
ylim([-1 0]);
xlim([0.125 2^8]);
set(gca,'XTick',dens,'XScale','log');
xlabel('Ripple density (cyc/oct)');
ylabel('Modulation depth');
axis square;

%% Reaction time / latency
figure(2)
for ii = 1:n
	subplot(3,ceil(n/3),ii)
	hist(L(ii,:),0:50:3000);
	xlim([300 2500]);
	xlabel('Reaction time (ms)');
	ylabel('Number of responses');
	str = ['Ripple density = ' num2str(dens(ii)) ' (cyc/oct)'];
	title(str)
	pa_verline(median(L(ii,:)));
end

%% Correct responses as a function of modulation depth
figure(3)
for ii = 1:n
	subplot(3,ceil(n/3),ii)
	plot(I(ii,:),R(ii,:)+0.01*randn(size(R(ii,:))),'k.','MarkerFaceColor','w');
	hold on
	axis square;
	ylim([-0.1 1.1]);
	pa_horline([0 1],'k:');
xlabel('Modulation depth');
ylabel('Response (+random jitter)');
	str = ['Ripple density = ' num2str(dens(ii)) ' (cyc/oct)'];
	title(str)

end

for ii = 1:n
	figure(3)
	subplot(3,ceil(n/3),ii)
	r = R(ii,:);
	i = I(ii,:);

	x = 0:2.5:100;
	nx = numel(x);
	
	NR = NaN(nx,1);
	for jj = 1:nx
		sel = i>=x(jj)-1.25 & i<=x(jj)+1.25;
		if sum(sel)
		NR(jj) = sum(r(sel))./sum(sel);
		end
	end
	sel = ~isnan(NR);
	
	plot(x(sel),NR(sel),'r-','LineWidth',2);
	hold on
	% plot(I,R,'ko','MarkerFaceColor','w');
	ylim([-0.1 1.1]);
	set(gca,'XScale','log','XTick',[1 10 100],'XTickLabel',[1 10 100]);
	axis square;
	% xlim([5 105]);
	xlim([1 100]);
	xlabel('Modulation depth (%)');
	ylabel('Probability Correct');
end
return
figure(1)
pa_datadir;
print('-dpng',mfilename);

