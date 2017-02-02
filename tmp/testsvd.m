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

nFreq   = 10;
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
M		= M+0.5*randn(size(M));
mu		= mean(M(:));
sd = std(M(:));
% mu = 0;
subplot(231);
imagesc(M')
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
	col = statcolor(64,[],[],[],'def',8,'disp',false);
colormap(col);

%% svd
[u,S,v] = svd((M-mu));
[n,p] = size(M);

pccoef = pca(M');

% keyboard
sigma = diag(S,0);

nsigma = numel(sigma);
Z = zeros(size(S));
rec = struct([]);
for ii = 1:nsigma
	Z(ii,ii) = sigma(ii);
	rec(ii).R =  u*Z*v'+mu;
	r = corrcoef(M(:),rec(ii).R(:));
	r = r(2)^2;
	rec(ii).r = r;

end


subplot(232)
imagesc(rec(1).R');
axis square;
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
r1 = corrcoef(M(:),rec(1).R(:));
title(num2str(r1(2)^2,'%0.2f'));


subplot(233)
imagesc(rec(2).R');
axis square;
axis square;
set(gca,'YDir','normal','TickDir','out');
box off
r2 = corrcoef(M(:),rec(2).R(:));
title(num2str(r2(2)^2,'%0.2f'));

%%
sigma = sigma.^2/(n-1);
alpha = sigma./sum(sigma);
% alpha	= sigma.^2/sum(sigma.^2);
beta	= pccoef/sum(pccoef);

subplot(234)
plot(alpha,'ko-','MarkerFaceColor','w');
hold on
plot(beta,'bs','MarkerSize',10);

for ii = 1:nsigma
	if ii>1
		r = rec(ii).r-rec(ii-1).r;
	else
		r = rec(ii).r;
	end
	plot(ii,r,'ro','MarkerSize',10);
end
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
for ii = 1:nsigma
	plot(ii,rec(ii).r,'ro','MarkerSize',10);
end

axis square;
box off
xlabel('Component #');
ylabel({'cumulative \Sigma ( \sigma^2/\Sigma(\sigma^2) )','R^2'});
ylim([-0.1 1.1]);
xlim([0 11]);
title('singular values \sigma');


subplot(236)
y = cumsum(alpha)-[rec.r]'

plot(y,'ko-','MarkerFaceColor','w');
hold on
% for ii = 1:nsigma
% 	plot(ii,rec(ii).r,'ro','MarkerSize',10);
% end
