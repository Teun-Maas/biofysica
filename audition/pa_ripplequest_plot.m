%% Initialization
close all
clear all
clc

%% Load the data
% pa_datadir('Ripplequest');
pa_datadir;
cd('Test');
fname = 'HH-HH-2012-02-02-0001-ripplequest';
fname = 'HH-IB-2012-02-02-0001-ripplequest';
fname = 'HH-MR-2012-02-02-0001-ripplequest';
fname = 'ripplequest';
fname = 'ripplequest_EGJ1';

% fname = 'ripplequest-MW-2012-02-04';
% fname = 'ripplequest-MW-2012-02-04-0001';
% fname = 'ripplequest-MW-2012-02-06-0001';

load(fname);
% Q(1).q
% return
%% Determine threshold based on responses
n		= numel(Q);
vel		= NaN(n,1);
dens	= NaN(n,1);
mod		= NaN(n,1);
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
	response		= q.response(1:ntrials);
	intensity		= 10.^(q.intensity(1:ntrials));
	I = [I intensity];
	R = [R response];
	L = [L lat];
	
	mod(ii)			= 10^t; % convert from logaritmic to linear
	vel(ii)			= Q(ii).vel;
	dens(ii)		= Q(ii).dens;
	fprintf('Final threshold estimate (mean±sd) is %.2f ± %.2f\n',mod,sd);


end

[dens,indx] = sort(dens);
mod = -mod(indx)/100;
vel = vel(indx);
plot(dens,mod,'ko-','MarkerFaceColor','w');
ylim([-1 0]);
set(gca,'XTick',dens,'XScale','log');
xlabel('Ripple density (cyc/oct)');
return
[I,indx] = sort(I);
R = R(indx);
L = L(indx);

% n = 1:numel(R);
% n = numel(R);
% C = cumsum(R)./n;
% C = C./max(C);

x = 0:2.5:100;
% x = pa_oct2bw(1,1:10);
sel = R==1;
Nall = hist(I,x);
N  = hist(I(sel),x);
N1 = N./Nall;
% N1 = 0.5*N./max(N);
N  = hist(I(~sel),x);
N0 = N./Nall;
% N0 = 0.5*N./sum(N);

N = [N1; N0];
isnan(N1)

XI	= I;
C	= interp1(x,N1,XI,'cubic');

figure
subplot(221)
% bar(x',N',1,'stack');
plot(I,C,'k-','LineWidth',2);
hold on
% plot(I,R,'ko','MarkerFaceColor','w');
ylim([-0.1 1.1]);
% set(gca,'XScale','log');
axis square;
% xlim([5 105]);
xlim([0 100]);
xlabel('Modulation depth (%)');
ylabel('Probability Correct');

subplot(222)
% plot(C,L,'k.');

sel = L<900;
X = C(sel);
Y = L(sel);
sigma = 0.1;
XI = min(C):0.01:max(C);
% keyboard
[mu,SE,XI] = pa_weightedmean(X,Y,sigma,XI);% set(gca,'XScale','log');
hold on
whos XI mu SE
pa_errorpatch(XI',mu',SE','r');
axis square
% xlim([1 99])
ylabel('Reaction time (ms)');
xlabel('Probability Correct');
drawnow
% return
%% Convert vectors to matrices
uvel		= unique(vel);
nvel		= numel(uvel);
udens		= unique(dens);
ndens		= numel(udens);
% sel = mod == min(mod);
% mod(sel) = NaN;

sel			= dens == 0 & vel==0;
mod(sel)	= nanmin(mod);
Z			= NaN(nvel,ndens);
for ii = 1:nvel
	for jj = 1:ndens
		sel = vel==uvel(ii) & dens == udens(jj);
		Z(ii,jj) = mod(sel);
	end
end
Z
% figure
subplot(223)
%% Threshold to sensitiy
Z = 100-Z;

% %% Graphics - simple truth
% h = imagesc(udens,uvel,Z);
% % image('YDir','reverse');
% get(h)
% % caxis([0 100]);
% axis square;
% set(gca,'Ydir','normal');
% colorbar
% caxis([40 100])
% % set(h);


% return
%% Interpolated/smoothed
X = udens;
Y = uvel;
XI = linspace(min(X),max(X),100);
YI = linspace(min(Y),max(Y),50)';
ZI = interp2(X,Y,Z,XI,YI,'linear');

%% Graphics - nice and cool
% figure
subplot(212)
contourf(XI,YI,ZI,0:100);
shading flat
h = colorbar;
%  = colorbar('peer',gca);
set(get(h,'ylabel'),'String', 'Sensitivity (%)');
% get(h)
% set(h,'YLabel','100-Modulation depth threshold (%)')
caxis([40 100]);
axis square;
xlabel('Ripple density (cycles/oct)');
ylabel('Ripple velocity (cycles/s)');
hold on
plot(0,1.8,'k^','MarkerFaceColor','k','MarkerSize',20);
colormap gray
% clabel('Modulation depth threshold (%)');
% pa_datadir;
% colormap('gray');
% print('-depsc','-painter',mfilename);

pa_datadir
% print('-depsc','MWtest');
print('-dpng','MWtest');