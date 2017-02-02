function avg = movcorr(data1,data2,windowSize)
% AVG = MOVCORR(DATA1,DATA2,WINDOWSIZE);
%
% Determines a moving average, AVG, of continuous DATA1 and DATA2
% with a window of WINDOWSIZE samples (default = 20).
%
% See also FILTE, MOVAVG
%
% https://nl.mathworks.com/matlabcentral/newsreader/view_thread/339632



%% Initialization
if nargin<1
	% just some test
    data    = normpdf(-100:100,30,30);
    noise1   = randn(1,length(data))/1000;
    noise2   = randn(1,length(data))/1000;
     data1    = [data+noise1 data+noise2];
   data2    = [data+noise2 -data-noise1];
end
if nargin<2
    windowSize = 20;
end

%%
x		= zscore(data1);
y		= zscore(data2);

x2		= x.^2;
y2		= y.^2;
xy		= x .* y;
A		= 1;
B		= ones(1,windowSize);
Stdx	= sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/windowSize))/(windowSize-1));
Stdy	= sqrt((filter(B,A,y2) - (filter(B,A,y).^2)*(1/windowSize))/(windowSize-1));
avg		= (filter(B,A,xy) - filter(B,A,x).*filter(B,A,y)/windowSize)./((windowSize-1)*Stdx.*Stdy);

avg = real(avg);
%% Some graphics if there is no other output
if nargout<1
    close all;
	subplot(211)
 h1 = plot(data1,'ko'); set(h1,'MarkerFaceColor','w');
    hold on;
 h2 = plot(data2,'ro'); set(h2,'MarkerFaceColor','w');
box off

 	subplot(212)
 h3 = plot(avg,'ko'); set(h3,'MarkerFaceColor','w');
ylim([-1.1 1.1]);
horline(0);
box off;
ylabel('Correlation');
%     h2 = plot(xavg,avg,'r-'); set(h2,'LineWidth',3);
%     xlabel('Sample Number');
%     ylabel('Amplitude');
end;