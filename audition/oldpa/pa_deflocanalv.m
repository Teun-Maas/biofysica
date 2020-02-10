function pa_deflocana(fname)

%% Clean
close all hidden
% clear all hidden


%% Load

if nargin<1
% 	cd('/Users/marcw/DATA/Student/Bahram/BY-SH-2014-10-03');
% 	fname = 'BY-SH-2014-10-03-0001';
% 	
% 	cd('/Users/marcw/DATA/tmp/MW-AR-2014-10-23');
% 	fname = 'MW-AR-2014-10-23-0001';
	cd('/Users/marcw/DATA/tmp/MW-PB-2014-11-19');
	fname = 'MW-PB-2014-11-19-0001';
	
	cd('/Users/marcw/DATA/tmp/default2snds/1snd');
	fname = 'RG-LV-2014-12-04-0001';
	
	cd('/Users/marcw/DATA/tmp/MW-RG-2014-12-05');
	fname = 'MW-RG-2014-12-05-0001';

% 	cd('/Users/marcw/DATA/tmp/LR-RG-2014-11-06');
% 	fname = 'LR-RG-2014-11-06-0001';

cd('/Users/marcw/DATA/tmp/RG-LV-2014-12-04');
fname = 'RG-LV-2014-12-04-0001';

cd('/Users/marcw/DATA/tmp/RG-LV-2014-12-16');
fname = 'RG-LV-2014-12-16-0001';

cd('/Users/marcw/DATA/tmp/lv');
fname = 'lv-0001';
else
	[pathstr,fname] = fileparts(fname);
	cd(pathstr);
end
load(fname) % load data
SupSac = pa_supersac(Sac,Stim,3,1); % combine saccade and stimulus matrix for first (1) sound (2) of trial
sel = ismember(SupSac(:,29),[65]);
SupSac = SupSac(sel,:);
%% Subdivide in sound types
sel		= ismember(SupSac(:,30),100:199); % BB sounds
BB		= SupSac(sel,:);
sel		= ismember(SupSac(:,30),200:299); % BB sounds
HP		= SupSac(sel,:);
sel		= ismember(SupSac(:,30),300:399); % BB sounds
LP		= SupSac(sel,:);


%% Graphics
% graphics(BB,0,1);
title({'';'Broad-band'}); % Provide all necessary information
% graphics(BB,1,4);
% graphics(HP,0,2);
title({'Azimuth';'High-pass'}); % Provide all necessary information
% graphics(HP,1,5);title({'Elevation';''});
graphics(LP,0,3);
title({'';'Low-pass'}); % Provide all necessary information
graphics(LP,1,6);


pa_datadir;
print('-depsc',mfilename);

function graphics(sac,dim,indx)

x	= sac(:,23+dim);
y	= sac(:,8+dim);
figure(1)
subplot(2,3,indx)
plot(x,y,'ko','MarkerFaceColor',[.7 .7 .7]);
axis square;							% publication-quality
box off;								% publication-quality
set(gca,'TickDir','out','TickLength',[0.005 0.025],...
	'XTick',-90:30:90,'YTick',-90:30:90,...
	'FontSize',10,'FontAngle','Italic','FontName','Helvetica'); % publication-quality
xlabel('Stimulus (deg)','FontSize',12);		% always provide an label on the x-axis
ylabel('Response (deg)','FontSize',12);	% and on the y-axis
axis([-90 90 -90 90]);						% set axis limits, this can be maximally between -90 and +90 deg
pa_unityline;

% Regression analysis
b = regstats(y,x,'linear',{'beta','rsquare','fstat'});
h = pa_regline(b.beta,'k-');
set(h,'linewidth',2);
b.fstat
if sign(b.beta(1)) >0
	str = ['\epsilon_R = ' num2str(round(b.beta(2)*10)/10) '\epsilon_T+' num2str(abs(round(b.beta(1)))) ];
else
	str = ['\epsilon_R = ' num2str(round(b.beta(2)*10)/10) '\epsilon_T-' num2str(abs(round(b.beta(1)))) ];
end

t = pa_text(0.1,0.9,str);
set(t,'FontSize',10)
