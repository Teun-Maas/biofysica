function double_rt

%% Initialization
close all

delay	= [0 5 10 20 40 80 120 160 320]';
ndelays = numel(delay);
Y = [];

for ii = 1:ndelays
	dT		= delay(ii);
	[y,x,RT] = getdata(1:9,dT);
			z		= (y-x(:,2))./(x(:,1)-x(:,2));
% 	a = (x(:,1)+x(:,2))/2;
% 	b  = (x(:,1)-x(:,2))/2;
% 	z		= ( y - a) ./ b ;

% 	RT = round(RT/25)*25;
% 	[uRT,~,sub] = unique(RT);
% 	
% 	mu = accumarray(sub,z,[],@mean);
% 	mu
	b = regstats(z,RT,'linear',{'beta','rsquare'});
	% 		figure(jj)
	figure(1)
	subplot(3,3,ii)
% 	plot(RT,z,'ko','MarkerFaceColor','w');
	hold on
% 	plot(uRT,mu)
	bubbleplot(RT,z,1,5);
	hold on
	x = [min(RT) max(RT)];
	y = b.beta(1)+x*b.beta(2);
	h = plot(x,y,'k-');
	
	set(h,'LineWidth',2);
	axis square;
	% 		ylim([-0.2 1.2]);
	xlim([300 800])
	horline([0 0.5 1],'k:');
	verline([0 delay(ii)]);
	text(700,0.9,num2str(b.rsquare,'%.2f'));
	ylim([-1 2]);
	
		x = 100:800;
	y = b.beta(1)+x*b.beta(2);

	Y = [Y;y];
	

end
savegraph('rt_delays','png');

% 	savegraph([mfilename num2str(jj)],'eps');
%%
col = statcolor(64,[],[],[],'def',8);
figure;
imagesc(Y)
caxis([0 1]);
colorbar
set(gca,'YDir','normal',...
	'TickDir','out',...
	'YTick',1:9,'YTickLabel',delay,...
	'XTIck',100:100:600,'XTickLabel',(100:100:600)+100)
colormap(col)
xlabel('Reaction time (ms)');
ylabel('Delay (ms)');
axis square
savegraph('rt_imagesc','png');
%%
% keyboard
return

%% The data specification:

delay	= [0 5 10 20 40 80 120 160 320]';
ndelays = numel(delay);
col = lines(9);
for jj = 1:9
	for ii = 1:ndelays
		dT		= delay(ii);
		[y,x,RT] = getdata(jj,dT);
		z		= (y-x(:,2))./(x(:,1)-x(:,2));
		b = regstats(z,RT,'linear',{'beta','rsquare'});
		% 		figure(jj)
		subplot(3,3,ii)
		plot(RT,z,'ko','MarkerFaceColor',col(jj,:));
		
		hold on
		x = [min(RT) max(RT)];
		y = b.beta(1)+x*b.beta(2);
		h = plot(x,y,'k-');
		
		set(h,'LineWidth',2)
		axis square;
		ylim([-0.2 1.2]);
		xlim([300 800])
		horline([0 0.5 1],'k:');
		verline([0 delay(ii)]);
		text(700,0.9,num2str(b.rsquare,'%.2f'));
	end
	% 	savegraph([mfilename num2str(jj)],'eps');
end


function [y,x,RT] = getdata(subject,dT)
load('/Users/marcw/Dropbox/Work/manuscripts/Rachel PhD/Papers/Paper 2 A9 Double Sounds/Data/A9DataMatrix.mat');


% 1) Trialnumber
% 2) 1snd location (el)
% 3) 2snd location (el)
% 4) el response
% 5) az resp (az location is always zero so not listed in matrix)
% 6) 1.snd
% 7) 2.snd
% 8) subjID (1-9)
% 9) delay (1:17 (-320:320)),
% 10) 1=single sound trials 2= bzz vs gwn trials 3= bzz&gwn (same location)


sel		= ismember(Matrix(:,10),[1 3]) & ismember(Matrix(:,8),subject);
y1		= Matrix(sel,4);
RT1		= Matrix(sel,11);
x1		= Matrix(sel,2);
N1		= numel(y1);
s1		= Matrix(sel,8);
us = unique(s1);
ns = numel(us);
for ii = 1:ns
	sel = s1==us(ii);
	b1(ii)		= regstats(y1(sel),x1(sel),'linear','beta');
end
delay	= -[-320 -160 -120 -80 -40 -20 -10 -5 0 5 10 20 40 80 120 160 320]';


sel			= ismember(delay(Matrix(:,9)),[dT -dT]) & Matrix(:,10)==2;
y			= Matrix(sel,4);
dt			= delay(Matrix(sel,9));
s			= Matrix(sel,8);
RT			= Matrix(sel,11);
x1			= Matrix(sel,2);
x2			= Matrix(sel,3);
x			= [x1 x2];
sel			= dt==-dT;
x(sel,:)	= [x2(sel) x1(sel)];

for ii = 1:ns
	sel = s==us(ii);
	sum(sel)
	x(sel,1) = b1(ii).beta(2)*x(sel,1)+b1(ii).beta(1);
	x(sel,2) = b1(ii).beta(2)*x(sel,2)+b1(ii).beta(1);
end

sel			= dt>=0;
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);
sum(sel)

sel			= ~isnan(RT);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

sel			= ismember(s,subject);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

d			= abs(x(:,2)-x(:,1));
sel			= d>20;
x			= x(sel,:);
y			= y(sel);
RT			= RT(sel);

isnan(RT)
hdi			= hdimcmc(RT);
% hdi = [200 700]
sel			= RT>hdi(1) & RT<hdi(2);
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
RT			= RT(sel);

N			= numel(y);


% b.beta
% x(:,1) = b.beta(2)*x(:,1)+b.beta(1);
% x(:,2) = b.beta(2)*x(:,2)+b.beta(1);



