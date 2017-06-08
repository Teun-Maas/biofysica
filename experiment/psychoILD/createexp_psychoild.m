function createexp_psychoild(varargin)
%
%
%
%

%% Initialization
x = 0:0.1:1;
x(x<0.001) = 0.001;
x(x>0.99) = 0.999;
niceild = probit(x);
niceild = niceild./max(abs(niceild));
ild		= keyval('ild',varargin,niceild); % ILDs

% Check this
nreps	= keyval('nreps',varargin,1); % number of repeats. NOTE: 1 = 0 repeats
bfreq	= keyval('bfreq',varargin,[750 2000 4000 6000]); % base frequency
shft	= keyval('shift',varargin, [-1 -1/3 -1/6 0 1/6 1/3 2/3 1]);
dispFlag	= keyval('display',varargin, true);


%% Combinatorics
[ild,bfreq,shft] = meshgrid(ild,bfreq,shft); % meshgrid

ild         = ild(:);

figure(2)
clf
x = -2:0.01:2;
[f,xi] = ksdensity(ild,x);
plot(xi,f)
hold on
N = hist(ild,x);
N = N./sum(N)*100;
plot(x,N)

freq_left   = bfreq(:);
freq_right  = oct2bw(bfreq(:),shft(:));

%% Repetition
ild         = repmat(ild,nreps,1);
freq_left	= repmat(freq_left,nreps,1);
freq_right	= repmat(freq_right,nreps,1);


%% Pseudo-randomize
n           = numel(ild);
idx         = randperm(n);
ild         = ild(idx);
freq_left	= freq_left(idx);
freq_right	= freq_right(idx);

ild = [ild; ild];
freq = [freq_left; freq_right];
freq_right = [freq_right; freq_left];
freq_left = freq;

dfreq       = freq_left-freq_right;
shft        = log2(freq_right./freq_left);

%% Get optimal ILD for frequency
% Check this!

i1 = 2*freq_left/1000;
i2 = 2*freq_right/1000;
mx = max([i1 i2],[],2);
ild = ild.*10;


%%
n           = numel(ild);

% cd('C:\DATA\SS\EXP');
save('psychoild1.exp','ild','freq_left','freq_right');

%% Display
if dispFlag
    figure(1)
    clf
    subplot(221)
    plot(freq_left,freq_right,'ko')
    %     set(gca,'XScale','log','YScale','log');
    set(gca,'XTick',unique(round(freq_left)),...
        'YTick',unique(round(freq_right)));
    axis([100 8000 100 8000]);
    axis square
    box off
    set(gca,'TickDir','out');
    xlabel('Frequency Left (Hz)');
    ylabel('Frequency Right (Hz)');
    title(n);
    
    subplot(222)
    plot(freq_left,ild,'ko')
    axis square
    box off
    set(gca,'TickDir','out');
    ylabel('ILD (dB)');
    xlabel('Frequency ');
    ylim([-50 50]);
    xlim([100 8000])
    
           subplot(223)
    plot(freq_left,shft,'ko')
    axis square
    box off
    set(gca,'TickDir','out');
    xlabel('frequency Left (Hz)');
    ylabel('Shift');
    xlim([100 8000]);
    ylim([-2.5 2.5])
    
        subplot(224)
    plot(ild,shft,'ko')
    axis square
    box off
    set(gca,'TickDir','out');
    xlabel('ILD (dB)');
    ylabel('Shift');
%     xlim([-30 30]);
    ylim([-2.5 2.5])
end