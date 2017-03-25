function pa_tonotopic
close all
clear all
clc;

f32 = 'E:\DATA\Cortex\Tonotopy\JoeData\';
xcl = 'E:\DATA\Cortex\Tonotopy';
joe(f32,xcl);


f32 = 'E:\DATA\Cortex\Tonotopy\ThorData\';
f32rm = 'E:\DATA\Cortex\Tonotopy\Data\';
xcl = 'E:\DATA\Cortex\Tonotopy\';
thor(f32,f32rm,xcl);

function joe(f32,xcl)
cd(xcl);
[NUMERIC,TXT] = xlsread('Joe-data.xlsx');
TXT		= TXT(2:end,:);
sel		= strcmpi(TXT(:,3),'tonerough');
NUMERIC = NUMERIC(2:end,:);
x		= NUMERIC(sel,5);
y		= NUMERIC(sel,6);
z		= NUMERIC(sel,1);
cd(f32)
files	= TXT(sel,1);
colormap jet
ls	= log10(logspace(log10(200),log10(12000),64));
k	= 0;
M	= NaN(length(files),4);
T	= NaN(length(files),1);
L	= NaN(length(files),1);
BW	= NaN(length(files),1);

for ii =1:length(files)
	fname = files{ii};
	if length(fname)>13
		DatFile = fname(1:6);
	elseif length(fname)==13
		DatFile= fname(1:5);
	end
	dname = [f32 DatFile,'experiment'];
	cd(dname);
	fname = pa_fcheckext(fname,'.f32');
	if exist(fname,'file')
			ToneChar = pa_tonechar(fname,'display',0);
			df = abs(ls-log10(ToneChar.charFrequency));

			k		= k+1;
			M(k,1)	= x(ii);
			M(k,2)	= y(ii);
			M(k,3)	= z(ii);
			M(k,4)	= ToneChar.charFrequency;
			
			if ~isempty(ToneChar.threshold)
			T(k) = ToneChar.threshold; % threshold
			L(k) = ToneChar.onsetLatency; % latency
			BW(k) = ToneChar.BW10; % bandwidth
			end
	else
		disp('Uh-oh');
		disp(fname);
	end
end

%%
sel		= isnan(M(:,1)) | isnan(M(:,2)) | isnan(M(:,4));
M		= M(~sel,:);
L		= L(~sel);
T		= T(~sel);
BW		= BW(~sel);

plotvor(M);
thresh	= [nanmean(T) nanstd(T)];
lat		= [nanmean(L) nanstd(L)];
bw		= [nanmean(BW) nanstd(BW)];
F = M(:,4);

figure(104)
subplot(131)
semilogx(F,T,'ko');
axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',0:25:125,'YTickLabel',0:25:125)
ylim([0 125]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('Threshold (dB SPL)');
str = ['Threshold = ' num2str(thresh(1),3) ' \pm ' num2str(thresh(2),2) ' (dB SPL)'];
title(str)

subplot(132)
semilogx(F,L,'ko');
axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',0:25:100,'YTickLabel',0:25:100)
ylim([0 100]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('Minimum latency (ms)');
str = ['Min Latency = ' num2str(lat(1),3) ' \pm ' num2str(lat(2),2)];
title(str)


subplot(133)
loglog(F,BW,'ko');
hold on
loglog([0.01 40]*1000,[3 3],'k:');

axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',[0.01 0.1 1 10],'YTickLabel',[0.01 0.1 1 10])
ylim([0.01 10]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('BW10dB');
str = ['BW10dB = ' num2str(bw(1),2) ' \pm ' num2str(bw(2),2)];
title(str)


%%
figure(4)
clf
CF			= M(:,4);
[uCF,~,IB]	= unique(CF);
nCF			= numel(uCF);
col			= jet(nCF);
col			= col(IB,:);
r1			= rand(size(M(:,1)))/2-0.25;
r2			= rand(size(M(:,2)))/2-0.25;

scatter(M(:,1)+r1,M(:,2)+r2,100,col,'filled')
hold on
errorbar(M(:,1)+r1,M(:,2)+r2,BW/2,'k.');
colorbar
axis square
axis([-3 8 -5 5])

%%
r1	= pa_freq2bw(200,CF);
r1	= pa_zscore(r1);
% r1 = r1-mean(r1);
r1 = (r1-min(r1))./(max(r1)-min(r1));
r1 = (r1-0.5)*2;
mf	= range(M(:,1)+r1)/range(M(:,1));
mu	= mean(M(:,1));
x	= (M(:,1)+r1-mu*mf)/mf+mu;

r2	= pa_freq2bw(200,CF);
r2	= pa_zscore(r2)/2;
r2 = (r2-min(r2))./(max(r2)-min(r2));
r2 = -(r2-0.5)*2;
mf	= range(M(:,2)+r2)/range(M(:,2));
mu	= nanmean(M(:,2));
y	= (M(:,2)+r2-mu*mf)/mf+mu;
uy = unique(M(:,2));
ny = numel(uy);
ux = unique(M(:,1));
nx = numel(ux);

figure(5)
scatter(x,y,100,col,'filled')
colorbar
axis square
axis([-3 8 -5 5])

x	= 100*x;
y	= 100*y;
xy	= [x y];
uxy = unique(xy,'rows');
nx	= length(uxy);
ucf = NaN(nx,1);
for ii = 1:nx
	sel = xy(:,1) == uxy(ii,1) & xy(:,2) == uxy(ii,2);
	ucf(ii) = median(pa_freq2bw(200,M(sel,4)));
end
z	= ucf;
x	= uxy(:,1);
y	= uxy(:,2);
XI	= linspace(min(x),max(x),10);
YI	= linspace(min(y),max(y),10);
[XI,YI] = meshgrid(XI,YI);
F	= TriScatteredInterp(x,y,z);
ZI	= F(XI,YI);

%%
figure(6)
% col = pa_statcolor(64,[],[],[],'def',2)
zcol = jet(64);
zcol =[1 1 1; zcol];
colormap(zcol);
imagesc(XI(1,:)/100,YI(:,1)/100,ZI);
shading flat;
cax = [0.5 4];
caxis(cax);
h = colorbar;
set(h,'YTick',1:cax(2),'YTickLabel',pa_oct2bw(200,(1:cax(2))));
axis square;
xlabel('Anterior-posterior');
ylabel('Medial-lateral');
axis([-3 8 -5 5])
box off
set(gca,'TickDir','out','YDir','normal');
pa_datadir
print('-depsc','-painter',[mfilename 'joe']);
%%
figure(7)
% col = pa_statcolor(64,[],[],[],'def',2)
zcol = jet(64);
zcol =[1 1 1; zcol];
% colormap(zcol);
% imagesc(XI(1,:)/100,YI(:,1)/100,ZI);
% shading flat;
% cax = [0.5 4];
% caxis(cax);
% h = colorbar;
% set(h,'YTick',1:cax(2),'YTickLabel',pa_oct2bw(200,(1:cax(2))));
% axis square;
% xlabel('Anterior-posterior');
% ylabel('Medial-lateral');
% axis([-3 8 -5 5])
% box off
% set(gca,'TickDir','out','YDir','normal');
% pa_datadir
% print('-depsc','-painter',[mfilename 'joe']);

% keyboard

%%
close all
% x = gallery('uniformdata',[10 2],5);
z	= ucf;
x	= uxy(:,1)/100;
y	= uxy(:,2)/100;
% x = [x y];
[~,~,zi] = unique(z);
col = jet(max(zi));
[v,c]=voronoin([x y]);
for i = 1:length(c)
	if all(c{i}~=1)   % If at least one of the indices is 1,
		% then it is an open region and we can't
		% patch that.

		patch(v(c{i},1),v(c{i},2),col(zi(i),:)); % use color i.
	end
end
axis([min(x) max(x) min(y) max(y)]);
%%
%%
% keyboard

function thor(f32,f32rm,xcl)
cd(xcl);
[NUMERIC,TXT]=xlsread('thor-data.xlsx');
TXT		= TXT(2:end,:);
sel		= strcmpi(TXT(:,3),'tonerough');
NUMERIC = NUMERIC(2:end,:);
x		= NUMERIC(sel,5);
y		= NUMERIC(sel,6);
z		= NUMERIC(sel,1);
cd(f32)
files	= TXT(sel,1);
colormap jet
c	= colormap;
ls	= log10(logspace(log10(200),log10(12000),64));
k = 0;
M = NaN(length(files),4);
T = NaN(length(files),1);
L = T;
BW = T;
for ii = 1:length(files)
	fname = files{ii};
	if length(fname)<15
		cd(f32)
		cd(fname(1:7))
	else
		cd(f32rm)
		cd(fname(1:end-5))
	end

	
	fname = pa_fcheckext(fname,'.f32');
	
	if exist(fname,'file')
		if ii>1
			figure(101)
			clf
		end
		
% 		try
		ToneChar = pa_tonechar(fname,'display',1);
% 								pause
		
if ~isnan(ToneChar.charFrequency)
		df = abs(ls-log10(ToneChar.charFrequency));
		[~,indx] = min(df);
		
		figure(666)
		plot(x(ii)+0.4*rand(1),y(ii)+0.4*rand(1),'o','Color',c(indx,:),'MarkerFaceColor',c(indx,:),'MarkerSize',5);
		hold on;
end

		k = k+1;
		M(k,1) = x(ii);
		M(k,2) = y(ii);
		M(k,3) = z(ii);
		M(k,4) = ToneChar.charFrequency;

		if ~isempty(ToneChar.threshold)
			T(k) = ToneChar.threshold; % threshold
			L(k) = ToneChar.onsetLatency; % latency
			BW(k) = ToneChar.BW10; % bandwidth
		end
	else
		disp('Uh-oh');
		disp(fname);
	end
end

%%
sel		= isnan(M(:,1)) | isnan(M(:,4));
M		= M(~sel,:);
L		= L(~sel);
T		= T(~sel);
BW		= BW(~sel);

plotvor2(M);
thresh	= [nanmean(T) nanstd(T)];
lat		= [nanmean(L) nanstd(L)];
bw		= [nanmean(BW) nanstd(BW)];
F = M(:,4);

figure(104)
subplot(131)
semilogx(F,T,'ko');
axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',0:25:125,'YTickLabel',0:25:125)
ylim([0 125]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('Threshold (dB SPL)');
str = ['Threshold = ' num2str(thresh(1),3) ' \pm ' num2str(thresh(2),2) ' (dB SPL)'];
title(str)

subplot(132)
semilogx(F,L,'ko');
axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',0:25:100,'YTickLabel',0:25:100)
ylim([0 100]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('Minimum latency (ms)');
str = ['Min Latency = ' num2str(lat(1),3) ' \pm ' num2str(lat(2),2)];
title(str)


subplot(133)
semilogx(F,BW,'ko');
hold on
semilogx([0.01 40]*1000,[3 3],'k:');
axis square;
set(gca,'XTick',[0.1 1 10 40]*1000,'XTickLabel',[0.1 1 10 40])
set(gca,'YTick',[0.01 0.1 1 10],'YTickLabel',[0.01 0.1 1 10])
ylim([0.01 10]);
xlim([0.1 40]*1000);
box off
xlabel('Characteristic frequency (kHz)');
ylabel('BW10dB');
str = ['BW10dB = ' num2str(bw(1),2) ' \pm ' num2str(bw(2),2)];
title(str)

lat = [nanmean(L) nanstd(L)]

%%
figure(3)
clf
CF = M(:,4);
[uCF,IA,IB] = unique(CF);
nCF = numel(uCF);
col = pa_statcolor(nCF,[],[],[],'def',2);
col = jet(nCF);
col = col(IB,:);

r1 = rand(size(M(:,1)))/2-0.25;
r2 = rand(size(M(:,2)))/2-0.25;

scatter(M(:,1)+r1,M(:,2)+r2,100,col,'filled')
hold on
errorbar(M(:,1)+r1,M(:,2)+r2,BW/2,'k.');

colorbar
axis square
axis([-3 8 -5 5])

%%
%%
r1	= pa_freq2bw(200,CF);
r1	= pa_zscore(r1);
% r1 = r1-mean(r1);
r1 = (r1-min(r1))./(max(r1)-min(r1));
r1 = (r1-0.5);
mf	= range(M(:,1)+r1)/range(M(:,1));
mu	= mean(M(:,1));
x	= (M(:,1)+r1-mu*mf)/mf+mu;

r2	= pa_freq2bw(200,CF);
r2	= pa_zscore(r2)/2;
r2 = (r2-min(r2))./(max(r2)-min(r2));
r2 = (r2-0.5)/2;
mf	= range(M(:,2)+r2)/range(M(:,2));
mu	= nanmean(M(:,2));
y	= (M(:,2)+r2-mu*mf)/mf+mu;
uy = unique(M(:,2));
ny = numel(uy);
ux = unique(M(:,1));
nx = numel(ux);

figure(7)
scatter(x,y,100,col,'filled')
colorbar
axis square
axis([-3 8 -5 5])

x	= 100*x;
y	= 100*y;
xy	= [x y];
uxy = unique(xy,'rows');
nx	= length(uxy);
ucf = NaN(nx,1);
for ii = 1:nx
	sel = xy(:,1) == uxy(ii,1) & xy(:,2) == uxy(ii,2);
	ucf(ii) = median(pa_freq2bw(200,M(sel,4)));
end
z	= ucf;
x	= uxy(:,1);
y	= uxy(:,2);
XI	= linspace(min(x),max(x),10);
YI	= linspace(min(y),max(y),10);
[XI,YI] = meshgrid(XI,YI);
F	= TriScatteredInterp(x,y,z);
ZI	= F(XI,YI);

figure(8)
% col = pa_statcolor(64,[],[],[],'def',2)
zcol = jet(64);
zcol =[1 1 1; zcol];
colormap(zcol);
imagesc(XI(1,:)/100,YI(:,1)/100,ZI);
shading flat;
cax = [0.5 4];
caxis(cax);
h = colorbar;
set(h,'YTick',1:cax(2),'YTickLabel',pa_oct2bw(200,(1:cax(2))));
axis square;
xlabel('Anterior-posterior');
ylabel('Medial-lateral');
axis([-3 8 -5 5])
box off
set(gca,'TickDir','out','YDir','normal');
pa_datadir
print('-depsc','-painter',[mfilename 'thor']);

%%
keyboard
function plotvor(M)
%%
sel = ~isnan(M(:,1)) & ~isnan(M(:,2)) & ~isnan(M(:,4));
M	= M(sel,:);

x	= M(:,1);
y	= M(:,2);
bf	= M(:,4);

xy		= [x y];
uxy		= unique(xy,'rows');
n		= length(uxy);
muBF	= NaN(n,1);
N = muBF;
colormap jet
c	= colormap;
ls	= log10(logspace(log10(200),log10(12000),64));
for ii = 1:n
	sel = x == uxy(ii,1) & y == uxy(ii,2);
	muBF(ii)	= nanmean(bf(sel));
	N(ii) = sum(sel);
	
	tmp = sort(bf(sel));
	ntmp = sum(sel);
	for jj = 1:ntmp
	figure(102)
	df = abs(ls-log10(tmp(jj)));
		[~,indx] = min(df);
	plot(uxy(ii,1),uxy(ii,1),'o','Color',c(indx,:),'MarkerFaceColor',c(indx,:),'MarkerSize',5);
	hold on
	end
	
end
axis square;
xlabel('Anterior-posterior');
ylabel('Medial-lateral');

y		= uxy(:,2);
x		= uxy(:,1);
y2		= (-6:0.5:-1)';
x2		= repmat(1.5,size(y2));
bf2		= repmat(100,size(y2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy = [x y];

y		= uxy(:,2);
x		= uxy(:,1);
y2		= (-6:0.5:-1)';
x2		= repmat(9,size(y2));
bf2		= repmat(100,size(y2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];

y		= uxy(:,2);
x		= uxy(:,1);
x2		= (2:0.5:9)';
y2		= repmat(-6.5,size(x2));
bf2		= repmat(100,size(x2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];


y		= uxy(:,2);
x		= uxy(:,1);
x2		= (2:0.5:9)';
y2		= repmat(-0.5,size(x2));
bf2		= repmat(100,size(x2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];

[uxy,indx]	= sortrows(uxy,[1,2]);
muBF		= muBF(indx);
[v,c] = voronoin(uxy);
figure(103)
subplot(121)
for j = 1:length(c)
	if all(c{j}~=1)   % If at least one of the indices is 1,then it is an open region and we can't patch that.
		hold on
		patch(v(c{j},1),v(c{j},2),log2(muBF(j))); % use color i.
	end
end
hold on
axis square;
cax = [8 13];
caxis(cax);
h = colorbar;
set(h,'YTick',cax(1):cax(2),'YTickLabel',2.^(cax(1):cax(2)));
axis square;
xlabel('Anterior-posterior');
ylabel('Medial-lateral');

xlim([2 8.5]);
xlabel('Anterior-posterior');
ylabel('Lateral-medial');
%%

function plotvor2(M)
%%
sel = ~isnan(M(:,1)) & ~isnan(M(:,2)) &~isnan(M(:,4));
M	= M(sel,:);

x	= M(:,1);
y	= M(:,2);
bf	= M(:,4);

xy		= [x y];
uxy		= unique(xy,'rows');
n		= length(uxy);
muBF	= NaN(n,1);
N		= muBF;
colormap jet
c	= colormap;
ls	= log10(logspace(log10(200),log10(12000),64));

for ii = 1:n
	sel = x == uxy(ii,1) & y == uxy(ii,2);
	muBF(ii)	= nanmean(bf(sel));
	N(ii) = sum(sel);
	
		tmp = sort(bf(sel));
	ntmp = sum(sel);
	for jj = 1:ntmp
	figure(102)
	df = abs(ls-log10(tmp(jj)));
		[~,indx] = min(df);
	plot(uxy(ii,1)+0.1*(jj-ntmp/2),uxy(ii,2),'o','Color',c(indx,:),'MarkerFaceColor',c(indx,:),'MarkerSize',5);
	hold on
	end
end


y		= uxy(:,2);
x		= uxy(:,1);
y2		= (-6:0.5:0)';
x2		= repmat(-3,size(y2));
bf2		= repmat(100,size(y2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy = [x y];

y		= uxy(:,2);
x		= uxy(:,1);
y2		= (-6:0.5:0)';
x2		= repmat(7,size(y2));
bf2		= repmat(100,size(y2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];

y		= uxy(:,2);
x		= uxy(:,1);
x2		= (-2:0.5:7)';
y2		= repmat(-6,size(x2));
bf2		= repmat(100,size(x2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];


y		= uxy(:,2);
x		= uxy(:,1);
x2		= (-2:0.5:7)';
y2		= repmat(2,size(x2));
bf2		= repmat(100,size(x2));
x		= [x; x2];
y		= [y; y2];
muBF	= [muBF; bf2];
uxy		= [x y];

[uxy,indx]	= sortrows(uxy,[1,2]);
muBF		= muBF(indx);
[v,c] = voronoin(uxy);
figure(103)
subplot(122)
for j = 1:length(c)
	if all(c{j}~=1)   % If at least one of the indices is 1,then it is an open region and we can't patch that.
		hold on
		patch(v(c{j},1),v(c{j},2),log2(muBF(j))); % use color i.
	end
end
hold on
% plot(uxy(:,1),uxy(:,2),'o');

axis square;
caxis([7 13])
h = colorbar;
set(h,'YTick',7:13,'YTickLabel',2.^(7:13));
xlabel('Anterior-posterior');
ylabel('Lateral-medial');
