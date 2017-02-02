
close all hidden
clearvars
clc

sdptf		= 1;
sdmtf		= 10;



%% PTF
vel		= 3;
dens	= 1;
mod		= 100; % Percentage (0-100%)
Fs		= 100;
mod		= mod/100; % Gain (0-1)
nTime	= 5;
time	= ((1:nTime)-1)/Fs; % Time (sec)

nFreq   = 11;
FreqNr  = 0:1:nFreq-1;
Oct     = FreqNr/20;                   % octaves above the ground frequency
oct		= repmat(Oct',1,nTime); % Octave

%% Generating the ripple
% Create amplitude modulations completely dynamic without static
A		= NaN(nTime,nFreq);
for ii = 1:nTime
	for jj = 1:nFreq
		A(ii,jj)      = 1 + mod*sin(2*pi*vel*time(ii) + 2*pi*dens*oct(jj)+0.5*pi);
	end
end

sd		= sdptf;
sdP		= sd/(2*pi);
sd		= sdmtf;
sdM		= sd/10;

M		= A;
M		= M-min(M(:));
M		= M./max(M(:));
M		= 2*pi*(M-0.5);
M		= M+sdP*randn(size(M));

M		= [flipud(M);M];
M = M+0.5*randn(size(M));
mu = mean(M(:));

% mu = 0;
subplot(231);
imagesc(M')
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
	col = statcolor(64,[],[],[],'def',8,'disp',false);
colormap(col);

%% svd
[u,S,v] = svd(M-mu);


% keyboard
sigma = diag(S,0);

S1		= zeros(size(S));
S1(1)	= sigma(1);
R1		= u*S1*v'+mu;

S2		= S1;
S2(2,2)	= sigma(2);
R2		= u*S2*v'+mu;

subplot(232)
imagesc(R1');
axis square;
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
r1 = corrcoef(M(:),R1(:));
title(num2str(r1(2)^2,'%0.2f'));


subplot(233)
imagesc(R2');
axis square;
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
r2 = corrcoef(M(:),R2(:));
title(num2str(r2(2)^2,'%0.2f'));

%%
alpha = sigma.^2/sum(sigma.^2);
subplot(234)
plot(alpha,'ko-','MarkerFaceColor','w');
hold on
plot(1,r1(2)^2,'ro','MarkerSize',10);
plot(2,r2(2)^2-r1(2)^2,'ro','MarkerSize',10);

axis square;
box off
xlabel('Component #');
ylabel('\sigma^2/\Sigma(\sigma^2) ');
ylim([-0.1 1.1]);
xlim([0 11]);
title('singular values \sigma');

subplot(235)
plot(cumsum(alpha),'ko-','MarkerFaceColor','w');
hold on
plot(1,r1(2)^2,'ro','MarkerSize',10);

plot(2,r2(2)^2,'ro','MarkerSize',10);
axis square;
box off
xlabel('Component #');
ylabel({'cumulative \Sigma ( \sigma^2/\Sigma(\sigma^2) )','R^2'});
ylim([-0.1 1.1]);
xlim([0 11]);
title('singular values \sigma');

%%
return
subplot(132)
imagesc(S);
axis square;

subplot(133)
imagesc(v);
axis square;

figure(3)
subplot(221)
plot(sigma,'ko-','MarkerFaceColor','w');
axis square;

