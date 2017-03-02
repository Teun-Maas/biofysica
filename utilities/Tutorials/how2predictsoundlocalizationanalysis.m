close all
clearvars
clc

X = -90:90;



%%
col		= lines(4);
ntar	= 15;
G		= [1 1.1 0.2 10];
STD		= [6 35 35 6];
ncond	= numel(G);
str = {'Normal-Hearing','Noisy','Noisy & Low Gain','Bilateral'};
xi = 0:100;
		n		= 100;

for Idx = 1:ncond
	for jj=1:4
		switch jj
			case 1
				r = 1;
				
				
			case 2
				r = 2;
			case 3
				r = 5;
				
			case 4
				r = 10;
				% 			sel = ismember(round(soundlevel),40:10:70)&ismember(soundtype,1:3);
		end
		% 	Tu = taraz(sel);
		Tu = 180*(rand(ntar*r,1)-0.5);
		Tu = round(Tu/15)*15;
		unique(Tu)
		g		= G(Idx);
		s		= STD(Idx);
		gain	= NaN(n,1);
		sd		= gain;
		mae		= gain;
		for ii = 1:n
			T = Tu;
			R = 180*(logisticfun(g*T,0,150,0.1)-0.5)+s*randn(size(T));
			b			= regstats(R,T,'linear',{'beta','r'});
			gain(ii)	= b.beta(2);
			sd(ii)		= std(b.r);
			N			= numel(T);
			mae(ii)		= mean(abs(R-T));
	
		end
		[MU,SD,A] = ellipse(sd,gain);


		subplot(ncond,6,6*(Idx-1)+1)
		plot(T,R,'k.');
		hold on;

		axis square;
		axis([-90 90 -90 90]);
		box off
		set(gca,'TickDir','out',...
			'XTick',-90:30:90,'YTick',-90:30:90);
		xlabel('Stimulus (deg)');
		ylabel('Response (deg)');
		title(str{Idx});
		unityline('k:');
		
		subplot(ncond,6,jj+6*(Idx-1)+1)
		plot(sd,gain,'k.','Color',col(jj,:));
		hold on
		plotellipse(MU,1.96*SD,A,'Color',col(jj,:));
		axis square;
		box off
		xlim([0 50]);
		ylim([-0.5 2]);
		horline(g,'k:');
		verline(s,'k:');
		title([num2str(N) ' trials']);
		plot(s,g,'ko','MarkerFaceColor','w','MarkerSize',10);
		xlabel('Standard deviation (deg)');
		ylabel('Gain');
		set(gca,'TickDir','out',...
			'XTick',0:10:50,'YTick',-1:0.5:2.5);
		
		
		subplot(ncond,6,6*(Idx-1)+6)

		f = ksdensity(mae,xi);

		plot(xi,f,'Color',col(jj,:));
		hold on
		axis square;
		xlim([0 80]);
		box off
		xlabel('Mean Absolute Error (deg)');
		ylabel('P');
				set(gca,'TickDir','out',...
			'XTick',0:10:50,'YTick',[]);
	end
end

savegraph(mfilename,'png');