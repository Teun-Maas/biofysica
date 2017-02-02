clearvars
close all
clc
cd('/Users/marcw/DATA/Peter Bremen');

load('ModData_hum.mat');


plot(Umd,mA');

%%
y = 1000./mA'/5;
% y = PC'/100;
x = round(100*Umd);
N = Wrong'+Correct';

% y = y(2:end,:);
% x = x(2:end);
% N = N(2:end,:);

y = floor(y.*N);
Y = [y(:) N(:)];
x = repmat(x,3,1);

s = repmat([1,2,3],length(y),1);
s = s(:);

sel = ismember(s,1:3);
psifit(x(sel),Y(sel,:),s(sel),'function',@logisticfun,'burnInSteps',10000,'gamma',0);

