%  psychometric function with large lapse rate
lambda = 0.3;
gamma = 0;
theta = -13;
omega = 10;
alpha = 0.1;
x = -30:10;
p = psifun(x,theta,omega,gamma,lambda,alpha);

% graphic visualization
close all
plot(x,p)
ylim([0 1]);

% get SNRS @ p = 25, 50 and 75% 
x = invpsifun([0.25 0.5 0.75],theta,omega,gamma,lambda,alpha)
% and check whether invpsifun really is the inverse (we should get 0.25
% 0.50 and 0.75]
p = psifun(x,theta,omega,gamma,lambda,alpha)

% and now determine SNRS @ 25, 50 and 75% of the 'curve'
x = invpsifun([0.25 0.5 0.75],theta,omega,0,0,alpha)
% and check at which p-values they are when considering guess and lapse
% rates
p = psifun(x,theta,omega,gamma,lambda,alpha)