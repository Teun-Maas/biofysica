close all

% x = gamrnd(1,1,300,4);
%
% y = skewness(x)
% xi = -5:0.5:5;
% for ii = 1:4
% 	N = hist(x(:,ii),xi)
% 	subplot(2,2,ii)
% plot(xi,N)
% title(y(ii))
% end
% return
fname = 'p14_celC.mat';


names = {
	'02'
	'04'
	'06'
	'08'
	'10'
	'14'
	'16'
	'18'
	'01'
	'03'
	'07'
	'09'
	'11'
	'13'
	'15'
	};


%%
lab = [zeros(1,8) ones(1,7)]+1;
labstr = {'C','B'};
freq = [60 50];
D = [];
% p = {'/Users/marcw/DATA/data pupil workshop 2015/Cel C/pupildata';'/Users/marcw/DATA/data pupil workshop 2015/Cel B/pupildata'};
xi = 2:0.03:9;
for i = 1:15
	path = which(['p',names{i},'_cel' labstr{lab(i)} '.mat']);
	load(path)
	figure(i)
	data	= pupilpreprocessing(rawdata,'freq',freq(lab(i)));
	D		= [D data];
	
	f		= ksdensity(data.smoothed(:),xi);
	f		= f./sum(f);
	
	figure(600)
	plot(xi,f+0.05*i);
	hold on
	x		= data.smoothed(:);
	x		= x(x>2);
	s		= skewness(x,0);
	[mx,x]	= max(f);
	text(xi(x),mx+i*0.05,num2str(s,'%.1f'),'HorizontalAlignment','center');
	text(2,i*0.05,names{i});
	% hist(data(:),xi);
	% axis([2 9 0 0.2]);
	% ylim([0 0.2])
	drawnow
end

%%

s = [D.skewness];
MPD = [D.MPD];
PPD = [D.PPD];

lat = [D.latency];

figure(700)
clf
subplot(121)
h = plot(s,MPD,'ko','MarkerFaceColor','w','MarkerSize',20);
axis square
box off
axis([-12 12 -0.08 0.28])
% lsline
axis([-12 12 -0.08 0.28])
set(gca,'TickDir','out','xtick',-10:5:10,'ytick',0:0.1:0.2);
horline;
verline;
xlabel('Skewness');
ylabel('Mean pupil dilation (mm)');
text(s,MPD,names,'HorizontalAlignment','center','VerticalAlignment','middle')
set(h,'MarkerEdgeColor','w');

% completely arbitrary "fit"
% assuming:
% 6 is strange outlier
% no/weak/random pupil effects when skewness/saturation is high
% stronger pupil effects when skewness/saturation is lower
x = -12:0.01:0;
p = 2.1*normpdf(x,0,4)+0.02;
hold on
plot(x,p,'k-')


subplot(122)
h = plot(s,PPD,'ko','MarkerFaceColor','w','MarkerSize',20);
axis square
box off
axis([-12 12 -0.08 0.68])
% lsline
axis([-12 12 -0.08 0.68])
set(gca,'TickDir','out','xtick',-10:5:10,'ytick',0:0.1:0.6);
horline;
verline;
xlabel('Skewness');
ylabel('Peak pupil dilation (mm)');
% text(s,PPD,names,'HorizontalAlignment','center','VerticalAlignment','middle')
text(s,PPD,labstr(lab),'HorizontalAlignment','center','VerticalAlignment','middle')
set(h,'MarkerEdgeColor','w');

x = -12:0.01:0;
p = 4*normpdf(x,0,4)+0.1;
hold on
plot(x,p,'k-')
% subplot(133)
% plot(s,lat,'ko','MarkerFaceColor','w')
% axis square
% box off
% % axis([-12 12 -0.08 0.28])
% lsline
% % axis([-12 12 -0.08 0.28])
% % set(gca,'TickDir','out','xtick',-10:5:10,'ytick',0:0.1:0.2);
% % horline;
% % verline;
% xlabel('Skewness');
% ylabel('Latency (s)');
% 

savegraph('pupilskewness','png');

%% SRT
fname		= which('SRT.xlsx');
[N,T,R]		= xlsread(fname);
SRT			= N(:,5);
xclnames	= N(:,1);
sel			= ~isnan(xclnames);
SRT			= SRT(sel);
xclnames	= xclnames(sel);

snames			= str2num(char(names));
[c,ia,ic]		= intersect(snames,xclnames);
x = SRT(ic);
y = PPD(ia)';
c = s(ia)';
sel = c>-10& y>0.06;
% sel = c>-5;
whos x y

b = regjags(y(sel),x(sel),'burnInSteps',2000,'numSavedSteps',2000);
b

	mubeta = [mean(b.beta0) mean(b.beta1)];

figure(900)
clf
subplot(121)
plot(x(sel),y(sel),'ko','MarkerFaceColor','w');
hold on
ylim([0 0.6]);
xlim([-20 5]);
col = [.7 .7 .7];
	for ii =  round(linspace(1,length(b.beta0),50))
		beta	= [b.beta0(ii) b.beta1(ii)]';
		h		= regline(beta,'k-');
		set(h,'Color',col);
	end
	
	h = regline(mubeta','k-');
set(h,'LineWidth',2);
plot(x(sel),y(sel),'ko','MarkerFaceColor','w','MarkerSize',15);
ylim([0 0.6]);
xlim([-20 5]);
axis square;
box off
set(gca,'TickDir','out');
xlabel('Speech reception threshold (dB)');
ylabel('Pupil dilation (mm)');

subplot(122)
plotpost(b.beta1(:));
axis square;