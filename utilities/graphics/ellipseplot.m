function h = pa_ellipseplot(Mu,SD,Phi,varargin)
% PA_ELLIPSEPLOT(MU,SD,A)
%
%  draw an ellipse with long and short axes SD(1) and SD(2)
%  with orientation A (in deg) at point Mu(1),Mu(2).
%
% see also PA_ELLIPSE

% 2011  Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
col         = keyval('Color',varargin);
if isempty(col)
    col = 'k';
end
Xo	= Mu(1);
Yo	= Mu(2);
L	= SD(1);
S	= SD(2);
DTR = pi/180;
Phi = Phi*DTR;

%% Ellipse
wt  = (0:.1:360).*DTR;
X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);

%% Graphics
h = patch(X,Y,col);
hold on
alpha(h,.2);
set(h,'EdgeColor',col,'LineWidth',2);

%% Rest
% wt = [0 180]*DTR;
% X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
% Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);
% plot(X,Y,'-','Color',Sty);
% 
% wt = [90 270]*DTR;
% X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
% Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);
% plot(X,Y,'-','Color',Sty);
