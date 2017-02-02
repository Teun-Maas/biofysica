function TOT = bubbleplot(X,Y,varargin)
% BUBBLEPLOT(X,Y)
%
% Make a bubbleplot of Y vs X
%
% BUBBLEPLOT(X,Y)

% (c) 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

%% Initialization

%%
[~,~,Xwidth,~] = binsize_scott(X);
[~,~,Ywidth,~] = binsize_scott(Y);
Xwidth = keyval('Xwidth',varargin,Xwidth);
Ywidth = keyval('Ywidth',varargin,Ywidth);


%% Histogram
X			= round(X/Xwidth)*Xwidth;
Y			= round(Y/Ywidth)*Ywidth;

% X		= round(X/afrnd)*afrnd;
% Y		= round(Y/afrnd)*afrnd;

uX			= unique(X);
uY			= unique(Y);

[UX,UY] = meshgrid(uX,uY);

%%
x		= uY;
TOT		= NaN(size(UX));
for i	= 1:length(uX)
	sel			= X == uX(i);
	r			= Y(sel);
	N			= hist(r,x);
	TOT(:,i)	= N;
end

%% Normalize
TOT		= log10(TOT+1);
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

SZ			= ceil(100*M);
[~,~,idx]	= unique(M);
% col			= flipud(gray(max(idx)));
col = statcolor(max(idx),[],[],[],'def',6);
% col			= parula(max(idx));
C			= col(idx,:);
h = scatter(x,y,SZ,C,'filled')

set(h,'MarkerEdgeColor','k');



