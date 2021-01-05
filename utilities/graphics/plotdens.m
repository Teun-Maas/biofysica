function [F,X1,X2] = plotdens(x,y,varargin)
% PLOTDENS(X,Y)
%
% Plot density
%
% See also KSDENSITY, CONTOURF

dep				= keyval('dependent',varargin,false); % scale y for each x
BW				= keyval('Bandwidth',varargin); % bandwidth for ksdensity
graph			= keyval('plot',varargin,true); % make plot
lim				= keyval('limits',varargin,[min(x) max(x) min(y) max(y)]); % make plot
ecol			= keyval('edgecolor',varargin,'k'); % scale y for each x
% cmap			= keyval('cmap'

X		= [x y];
% if dep
gridx1	= linspace(lim(1),lim(2),100);
gridx2	= linspace(lim(3),lim(4),100);

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
	if exist('cbrewer.m','file')
% 		cmap	= cbrewer('seq','Reds',64,'pchip');
		cmap	= flipud(cbrewer('div','RdBu',64,'pchip'));
		
	else
		cmap	= statcolor(64,[],[],[],'def',6);
	end
	[~,h] = contourf(X1,X2,F,0:0.01:1);
	h.EdgeColor = 'none';
	hold on
	[~,h] = contour(X1,X2,F,0:0.1:1);
	h.EdgeColor = ecol;
	colormap(cmap)
	% colorbar
	set(gca,'Ydir','normal');
	hold on
end