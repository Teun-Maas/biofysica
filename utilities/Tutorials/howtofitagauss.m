function pa_howtofitagauss
% PA_HOWTOFITAGAUSS
%
% Example/demo script for fitting a Gauss
%
% See also NLINFIT, PA_GAUSS

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Cleaning
close all
clear all
clc

%% Produce gaussian data
x = 0:.1:20;
b = [10 2 6 1]; % parameters mean, sd, amplitude and offset
y = pa_gauss(b,x); % make gaussian data
plot(x,y,'ko-')
hold on
y = y+randn(size(x)); % add some noise
% Graph
plot(x,y,'bo-')

%% Let's fit a gaussian through the data
[mx,indx] = max(y); % Estimate some parameters
beta		= nlinfit(x,y,@pa_gauss,[x(indx) 1 mx 0]); % Do a non-linear gaussian fit (save the function gauss below as a separate m-file)
% Check how well fitting was
y2			= pa_gauss(beta,x);
[r,p]		= corrcoef(y,y2);
r = r(2)^2;
p = p(2);
% plot the fitted gaussian
plot(x,y2,'r-');

title([r p]);