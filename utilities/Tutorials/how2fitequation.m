function tmp
%% 
close all
clearvars


%% Simulate some data

%% Stimuli
x	= -90:90; % locations
udx	= -60:10:60; % spatial disparity
uadx = unique(abs(udx));
nuadx = numel(uadx);
ndx = numel(udx);
udt = [0 5 10 20 40 80]; % temporal disparity
udt = 0:80;
ndt = numel(udt);
[x,dx,dt] = ndgrid(x,udx,udt);
x	= x(:);
dx	= dx(:);
dt = dt(:);

x1	= x; % location 1
x2	= x+dx; % location 2

sel = abs(x2)<90; % remove all non-existent locations
x2	= x2(sel);
x1	= x1(sel);
dx	= dx(sel);
dt = dt(sel);

%% Responses
g			= 1; % response gain is optimal
b			= 0; % response bias is zero
sd			= 8; % additive response noise

%  % linear dependence of spatial disparity
% g_dx		= 1; % how weights depend on spatial disparity - slope
% b_dx		= 0; % how weights depend on spatial disparity -offset
% omega		= g_dx*(90-abs(dx))+b_dx;


%% non-linear dependence of spatial disparity
g_dx		= 0.01; % how weights depend on spatial disparity - slope
b_dx		= 0; % how weights depend on spatial disparity -offset
omega		= g_dx*(90-abs(dx)).^2+b_dx;
% try to figure out what the omega-dx relationship should be 
% - linear is the simplest
% - quadratic might be nice, perhaps power, perhaps exponential
% whatever you decide to use, you also need to adapt the fitted equation
% function @fun

%%

a			= 0.05*2;
z			= 2*log(1/a-1); % in this parametrization, omega is the delay for which the weight is 0.95
w			= 1./(1+exp(-z./omega.*dt)); % weights depend on temporal disparity in a sigmoidal fashion

[~,~,subs]	= unique([abs(dx) dt],'rows');
W			= accumarray(subs,w,[],@mean); % average weight per dt-dx combination - without noise this is trivial
W			= reshape(W,ndt,numel(W)/ndt); % reshaping to a matrix for easy plotting

R			= g*(w.*x1 + (1-w).*x2) + b + sd.*randn(size(x1));

figure(1)
clf

subplot(231)
% plot(dt,w,'o','MarkerFaceColor','w');
hold on
h=plot(udt,W,'-','MarkerFaceColor','w');
nicegraph;
legend(h,num2str(uadx'),'Location','SE');

subplot(234)
plot(x1,R,'.')
xlabel('delay (ms)');
ylabel('weight');


subplot(235)
plot(x2,R,'.')

subplot(236)
plot((x1+x2)/2,R,'.')


for ii = 4:6
	subplot(2,3,ii)
	nicegraph
	xlabel('target (deg)');
	ylabel('response (deg)');
end

%% Regression
x		= [x1 x2 dt dx]; % independent variables
y		= R; % dependent variable
beta0	= [1 0 0.01 0]; % initial values
[beta,resid,~,covb]	= nlinfit(x,y,@fun,beta0);


%% Robust fitting - typically not needed
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare'; % 
% beta0 = [1 0 1 0];
% [beta,resid,~,covb]	= nlinfit(x,y,@fun,beta0,opts);

%% Confidence intervals on parameters to check whether you are confident on the fit
% ci			= nlparci(beta,resid,'covar',covb);

%% Prediction - to check whether the model correctly predicts/describes the data
% [ypred,delta] = nlpredci(@fun,x,beta,resid,'covar',covb);

%% fitted weights

g_dx		= beta(3); % how weights depend on spatial disparity - slope
b_dx		= beta(4); % how weights depend on spatial disparity -offset
figure(1)
subplot(232)
hold on
for ii = 1:nuadx
omega		= g_dx*(90-uadx(ii)).^2+b_dx;
a			= 0.05*2;
z			= 2*log(1/a-1); % in this parametrization, omega is the delay for which the weight is 0.95
w			= 1./(1+exp(-z./omega.*udt)); % weights depend on temporal disparity in a sigmoidal fashion
plot(udt,w);
end
nicegraph
legend(num2str(uadx'),'Location','SE');


%% Weight-dx
omega		= g_dx*(90-uadx).^2+b_dx;

subplot(233)
plot(uadx,omega,'o-','MarkerFaceColor','w');
nicegraph
xlabel('spatial disparity (deg)');
ylabel('omega');

%%


function y = fun(beta,x)

x1			= x(:,1);
x2			= x(:,2);
dx			= x(:,4);
dt			= x(:,3);
g			= beta(1); % response gain is optimal
b			= beta(2); % response bias is zero
g_dx		= beta(3); % how weights depend on spatial disparity - slope
b_dx		= beta(4); % how weights depend on spatial disparity -offset
omega		= g_dx*(90-abs(dx)).^2+b_dx;
a			= 0.05*2;
z			= 2*log(1/a-1); % in this parametrization, omega is the delay for which the weight is 0.95
w			= 1./(1+exp(-z./omega.*dt)); % weights depend on temporal disparity in a sigmoidal fashion
y			= g*(w.*x1 + (1-w).*x2) + b;



