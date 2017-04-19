
function tmp

% Yamada, T., Umeyama, S., Matsuda, K., Tanikawa, Y., & Yamada, Y. (2012)
% Separation of fNIRS Signals into Functional and Systemic Components Based
% on Differences in Hemodynamic Modalities. PLoS One, 7, e50271. 
kf = -0.6; % universal 

a = rand(100,1);
b = rand(100,1);
[A,B] = meshgrid(a,b);
B = 2*A+B;
M = [A(:) B(:)];

close all
plot(M(:,1),M(:,2),'.');

MI = mutualinformation(A(:),B(:))

title(MI)
%%

%%
% keyboard
function MI = mutualinformation(dHbOf,dHbOs)

f			= dHbOf;
s			= dHbOs;
uf			= unique(f);
us			= unique(s);
ufs			= unique([f s],'rows');
[pf,xf]		= ksdensity(f,uf);
[ps,xs]		= ksdensity(s,us);
[pfs,xfs]	= ksdensity([f s],ufs);

% pfs = ksdensity?

MI = 0;
for idx = 1:size(xfs,1)
	self = xf(:,1) == xfs(idx,1);
	sels = xs(:,1) == xfs(idx,2);
	
MI = MI+pfs(self).*log( pfs(idx)./(pf(self).*ps(sels)) );
end


function [F,S] = functionalvssystemic(dHbO,dHbr)
F = 1/(kf-ks)*[-ks 1; -kf*ks kf]*[dHbO;dHbr];
% dHbOf = F(1,:);
% dHbRf = F(2,:);
S = 1/(kf-ks)*[ks 1; kf*ks -ks]*[dHbO;dHbr];
% dHbOs = S(1,:);
% dHbRs = S(2,:);