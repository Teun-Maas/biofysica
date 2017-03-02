function pa_tmp

close all hidden
clear all hidden
clc

sdptf		= 1;


%% PTF
vel		= 4;
dens	= 0.5;
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

PTF		= A;
PTF		= PTF-min(PTF(:));
PTF		= PTF./max(PTF(:));
PTF		= 2*pi*(PTF-0.5);
PTF		= PTF+sdP*randn(size(PTF));


figure(666)
imagesc(PTF')
axis square;
col = pa_statcolor(64,[],[],[],'def',8,'disp',false);
colormap(col);
set(gca,'XTick',[],'YTick',[]);
xlabel('X');
ylabel('Y');
title('HCL Divergent');
% colorbar

print('-depsc','-painter',[mfilename 'hcl1']);

figure(667)
imagesc(PTF')
axis square;
set(gca,'XTick',[],'YTick',[]);

col = jet(64);
colormap(col);
xlabel('X');
ylabel('Y');
title('RGB');
print('-depsc','-painter',[mfilename 'jet']);

% for ii = 2:9
% figure(667+ii)
% imagesc(PTF')
% axis square;
% col = pa_statcolor(64,[],[],[],'def',ii,'disp',false);
% % col = flipud(col);
% colormap(col);
% set(gca,'XTick',[],'YTick',[]);
% xlabel('X');
% ylabel('Y');
% title('HCL Rainbow');
% print('-depsc','-painter',[mfilename 'hcl' num2str(ii)]);
% end
% colorbar

