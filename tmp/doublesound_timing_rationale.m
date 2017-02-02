% t = 0:10;
%
% y = [0.5 1 1 1 1 1

clearvars;
close all

%%
	subplot(1,2,1)
	plot([0 1],[0.5 1],'k-');
	hold on
	plot([1 4],[1 1],'k-');
	x = 4:0.1:4.9;
	y = +1-0.1*(x-4);
	plot(x,y,'k.');
	
	y = +1-0.9*(x-4);
	plot(x,y,'k.');
	
	plot([5 10],[1 1],'k-');
	
	plot([5 10],[0 0],'k-');
	
	ylim([-0.1 1.1]);
	
	xlim([-1 11]);
	
	axis square
	box off
	verline([0 1 5]);
	plot([4 4],[-1 0.2],'k--');
	plot([4 4],[0.8 1.1],'k--');
	
	text(0.5,0,'Time-Weighted Average',...
		'HorizontalAlignment','left',...
		'VerticalAlignment','middle',...
		'Rotation',90);
	
	text(2.5,0,'Precedence Effect',...
		'HorizontalAlignment','left',...
		'VerticalAlignment','middle',...
		'Rotation',90);
	
	text(4,0.5,'Echo Threshold',...
		'HorizontalAlignment','center',...
		'VerticalAlignment','middle',...
		'Rotation',90);
	
	text(4.5,0,'Bimodal Double TWA',...
		'HorizontalAlignment','left',...
		'VerticalAlignment','middle',...
		'Rotation',90);
	
	text(7.5,0.05,'Bimodal Double',...
		'HorizontalAlignment','left',...
		'VerticalAlignment','middle',...
		'Rotation',90);
	xlabel('Time (ms)');
	ylabel('Perceived location');
	set(gca,'YTick',[0 1],'YTickLabel',{'Lag','Lead'});
	title('Azimuth');

%%
	subplot(1,2,2)
hold on


ms  = 0.6+(0:10)*0.04;
ms = 20*ms;
ms2 = 21-ms;
% for ii = 0:10
% 	plot(ii*30,1,'ko','MarkerFaceColor','w','MarkerSize',ms2(ii+1));
% 	plot(ii*30,0,'ko','MarkerFaceColor','w','MarkerSize',ms2(ii+1));
% 	
% 	plot(ii*30,0.5+ii*0.05,'ko','MarkerFaceColor','w','MarkerSize',ms(ii+1));
% 
% end
% plot([0 300],[1 1],'k:');
% plot([0 300],[0 0],'k:');

plot([0 300],[0.5 1],'k-');
plot([300 350],[1 1],'k-');

axis square;
ylim([-0.1 1.1])
	xlabel('Time (ms)');
	ylabel('Perceived location');
	set(gca,'YTick',[0 1],'YTickLabel',{'Lag','Lead'});
	title('Elevation');
xlim([-10 350])
verline([0 300]);
	text(150,0.5,'Time-Weighted Average',...
		'HorizontalAlignment','center',...
		'VerticalAlignment','middle');
	
%%
savegraph('rationaleandconclusion','png');