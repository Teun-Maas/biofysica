%% Clear
close all hidden
clear all hidden
clc

%% Just some random data
x = rand(5,10);

%% 2 dimensions
% Plot first 2 dimensions in 2D space, X vs Y
subplot(141)
plot(x(1,:),x(2,:),'ko','MarkerFaceColor','w');
axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('2 dimensions');
xlabel('X');
ylabel('Y');
axis([-0. 1.1 -0.1 1.1]);

%% 3 dimensions
% Add 3rd dimension as colour
c = round(x(3,:)*10);
c = (c-min(c))+1;
nc = max(c);
n = numel(c);
if pa_isodd(nc)
	nc = nc+1;
end
col = pa_statcolor(nc,[],[],[],'def',8);
subplot(142)
for ii = 1:n
	indx = c(ii);
	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8]);
	hold on
end
axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('3 dimensions');
xlabel('X');
ylabel('Y');
axis([-0. 1.1 -0.1 1.1]);

%% 4 dimensions
% 4th dimension is coded as marker-size
subplot(143)
s = 2*round(x(4,:)*10);
s = (s-min(s))+1;
for ii = 1:n
	indx = c(ii);
	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
	hold on
end
axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('4 dimensions');
xlabel('X');
ylabel('Y');
axis([-0. 1.1 -0.1 1.1]);

%% 5 dimensions
% the 5th dimension is coded by symbol
subplot(144)
g = x(5,:)<0.5;
for ii = 1:n
	indx = c(ii);
	if g(ii)
	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
	else
	plot(x(1,ii),x(2,ii),'kd','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
	end
	hold on
end
axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('5 dimensions');
xlabel('X');
ylabel('Y');
axis([-0. 1.1 -0.1 1.1]);

%% Save
print('-depsc','-painter',mfilename);


%%
figure
%% 2 dimensions
% Plot first 2 dimensions in 2D space, X vs Y
subplot(141)
bar(x(1,:),'FaceColor',[.7 .7 .7]);
axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('2 dimensions');
xlabel('X');
ylabel('Y');
axis([0 size(x,2)+1 0 1.1*max(x(1,:))]);

subplot(142)
bar3(x,1);
% axis square;
box off
set(gca,'XTick',[],'YTick',[]);
title('3 dimensions');
xlabel('X');
ylabel('Y');
% axis([0 size(x,2)+1 0 1.1*max(x(1,:))]);

% %% 3 dimensions
% % Add 3rd dimension as colour
% c = round(x(3,:)*10);
% c = (c-min(c))+1;
% nc = max(c);
% n = numel(c);
% if pa_isodd(nc)
% 	nc = nc+1;
% end
% col = pa_statcolor(nc,[],[],[],'def',8);
% subplot(142)
% for ii = 1:n
% 	indx = c(ii);
% 	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8]);
% 	hold on
% end
% axis square;
% box off
% set(gca,'XTick',[],'YTick',[]);
% title('3 dimensions');
% xlabel('X');
% ylabel('Y');
% axis([-0. 1.1 -0.1 1.1]);
% 
% %% 4 dimensions
% % 4th dimension is coded as marker-size
% subplot(143)
% s = 2*round(x(4,:)*10);
% s = (s-min(s))+1;
% for ii = 1:n
% 	indx = c(ii);
% 	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
% 	hold on
% end
% axis square;
% box off
% set(gca,'XTick',[],'YTick',[]);
% title('4 dimensions');
% xlabel('X');
% ylabel('Y');
% axis([-0. 1.1 -0.1 1.1]);
% 
% %% 5 dimensions
% % the 5th dimension is coded by symbol
% subplot(144)
% g = x(5,:)<0.5;
% for ii = 1:n
% 	indx = c(ii);
% 	if g(ii)
% 	plot(x(1,ii),x(2,ii),'ko','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
% 	else
% 	plot(x(1,ii),x(2,ii),'kd','MarkerFaceColor',col(indx,:),'MarkerEdgeColor',[.8 .8 .8],'MarkerSize',s(ii));
% 	end
% 	hold on
% end
% axis square;
% box off
% set(gca,'XTick',[],'YTick',[]);
% title('5 dimensions');
% xlabel('X');
% ylabel('Y');
% axis([-0. 1.1 -0.1 1.1]);

%% Save
print('-depsc','-painter',mfilename);
