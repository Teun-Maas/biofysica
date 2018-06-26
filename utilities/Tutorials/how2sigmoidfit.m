function how2sigmoidfit


close all
clearvars
x		= -90:90;
G		= 1:0.1:5; % linear target-response gain
nG = numel(G);
D		= NaN(size(G)); % derivative at T=0 deg
C		= D; % compression
B		= D; % response bias (deg)
T		= D; % target bias
O		= D;
for ii = 1:nG
	g		= G(ii);
	y		= g*x+0*randn(size(x));
	y		= tanh(y/90)*90;
% 		y = y*0.7+10;
	beta0	= [0 1 1 0];
	beta	= nlinfit(x,y,@sigmoidfun,beta0);
	ypred	= sigmoidfun(beta,x);
	
	
	figure(1)
	clf
	subplot(121)
	plot(x,y,'.');
	hold on
	plot(x,ypred,'o');
	nicegraph
	unityline
	xlabel('Target (deg)');
	ylabel('Response (deg)');
	
	
	C(ii) = beta(3);
	B(ii) = beta(4);
	O(ii) = beta(2);
	T(ii) = beta(1);
	%% Derivative
	% http://math2.org/math/derivatives/more/hyperbolics.htm
	d		= beta(3)*beta(2)*(1-tanh(beta(2)*x/90).^2); % 'full' derivative at all T
	D(ii) = beta(3)*beta(2); % derivative at T=0
	
	
	figure(1)
	subplot(122)
	plot(x,d)
	nicegraph
	
	drawnow
end

%%

figure(2)
clf
plot(G,D)
nicegraph
xlabel('Linear target-response gain');
ylabel('Derivative at \alpha_T = 0 deg');

figure(3)
clf
subplot(221)
plot(G,T)
nicegraph
xlabel('Linear target-response gain');
ylabel('Target bias (deg)');

subplot(222)
plot(G,O)
nicegraph
xlabel('Linear target-response gain');
ylabel('\omega');

subplot(223)
plot(G,C)
nicegraph
xlabel('Linear target-response gain');
ylabel('Compression factor');

subplot(224)
plot(G,B)
nicegraph
xlabel('Linear target-response gain');
ylabel('Response bias deg');


function y = sigmoidfun(beta,x)
t		= beta(1); % target bias
g		= beta(2); % gain
c		= beta(3); % compression 
b		= beta(4); % response bias
y		= 90*c*tanh((g*(x-t))/90)+b;
