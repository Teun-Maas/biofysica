function densplot(x,y,varargin)
% DENSPLOT(X,Y)
%
% Plot density
%
% See also KSDENSITY, CONTOURF

dep				= keyval('dependent',varargin,false); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
BW				= keyval('Bandwidth',varargin); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
lim				= keyval('limit',varargin); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun

X		= [x y];
if isempty(lim)
	lim = [min(x) max(x) min(y) max(y)];
end
gridx1	= linspace(lim(1),lim(2),100);
gridx2	= linspace(lim(3),lim(4),100);
[x1,x2] = meshgrid(gridx1, gridx2);
[m,n]	= size(x1);
x1		= x1(:);
x2		= x2(:);
xi		= [x1 x2];
if isempty(BW)
	[F,XI] =  ksdensity(X,xi);
else
	[F,XI] =  ksdensity(X,xi,'Bandwidth',BW);
end
F		= reshape(F,m,n);
F		= F./sum(F(:));
if dep
	mu		= max(F);
	F = bsxfun(@rdivide,F,mu);
end

X1		= reshape(XI(:,1),m,n);
X2		= reshape(XI(:,2),m,n);
cmap	= statcolor(64,[],[],[],'def',6);
[~,h] = contourf(X1,X2,F.^2,0:0.01:1);
h.EdgeColor = 'none';
hold on
colormap(cmap)
% colorbar
set(gca,'Ydir','normal');
hold on