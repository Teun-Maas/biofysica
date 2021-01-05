function [TOT,s] = bubbleplot(X,Y,varargin)
% BUBBLEPLOT(X,Y)
%
% Make a bubbleplot of Y vs X.
%
% BUBBLEPLOT(X,Y)
%
% - Xwidth
% - Ywidth
% - col
% - dependent
% - MarkerSizeFactor

% (c) 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

%% Initialization

%%
[~,~,Xwidth,~] = binsize_scott(X);
[~,~,Ywidth,~] = binsize_scott(Y);
Xwidth			= keyval('Xwidth',varargin,Xwidth);
Ywidth			= keyval('Ywidth',varargin,Ywidth);

def				= keyval('col',varargin,6);
dep				= keyval('dependent',varargin,false); % scale y for each x
MSF				= keyval('MarkerSizeFactor',varargin,3); % scale y for each x

%% Histogram
X			= round(X/Xwidth)*Xwidth;
Y			= round(Y/Ywidth)*Ywidth;

% X		= round(X/afrnd)*afrnd;
% Y		= round(Y/afrnd)*afrnd;

uX			= unique(X);
uY			= unique(Y);

[UX,UY] = meshgrid(uX,uY);

% if numel(uY)==1
% 	keyboard
% end

%%
x		= uY;

TOT		= NaN(size(UX));
for ii	= 1:length(uX)
	sel			= X == uX(ii);
	r			= Y(sel);
	N			= hist(r,x);
	if isscalar(x)
		N			= histogram(r,[x-Xwidth x+Xwidth]);
	end
	TOT(:,ii)	= N;
end


if dep
	mu		= max(TOT);
	TOT = TOT./mu;
end

%% Normalize
% TOT		= log10(TOT+1);
mxTOT	= nanmax(nanmax(TOT));
mnTOT	= nanmin(nanmin(TOT));
TOT		= (TOT-mnTOT)./(mxTOT-mnTOT);

%% Plot
M	= TOT(:);
x	= UX(:);
y	= UY(:);

sel = M>0;
M	= M(sel);
x	= x(sel);
y	= y(sel);

SZ			= ceil(MSF*100*M);
[~,~,idx]	= unique(M);
try
	col = cbrewer('seq','OrRd',max(idx),'pchip');
catch
	col			= statcolor(max(idx),[],[],[],'def',def);
end
C			= col(idx,:);
s			= scatter(x,y,SZ,C,'filled');
set(s,'MarkerEdgeColor','k');



