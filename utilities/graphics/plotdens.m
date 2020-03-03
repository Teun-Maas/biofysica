function [F,X1,X2] = plotdens(x,y,varargin)
% PLOTDENS(X,Y)
%
% Plot density
%
% See also KSDENSITY, CONTOURF

dep				= keyval('dependent',varargin,false); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
BW				= keyval('Bandwidth',varargin); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
graph			= keyval('plot',varargin,true); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun

X		= [x y];
% if dep
gridx1	= linspace(min(x),max(x),100);
gridx2	= linspace(min(y),max(y),100);

% else
% gridx1	= linspace(prctile(x,5),prctile(x,95),100);
% end
% gridx2	= linspace(prctile(y,5),prctile(y,95),100);
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
F		= F./max(F(:));
if dep
	mu		= max(F);
	F = F./mu;
% 	F = bsxfun(@rdivide,F,mu);
end

X1		= reshape(XI(:,1),m,n);
X2		= reshape(XI(:,2),m,n);
if graph
	cmap	= statcolor(64,[],[],[],'def',6);
	[~,h] = contourf(X1,X2,F,0:0.01:1);
	h.EdgeColor = 'none';
	hold on
	colormap(cmap)
	% colorbar
	set(gca,'Ydir','normal');
	hold on
end