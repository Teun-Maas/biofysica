function [F,S] = functionalvssystemic(dHbO,dHbR,kf,ks)

% Yamada, T., Umeyama, S., Matsuda, K., Tanikawa, Y., & Yamada, Y. (2012)
% Separation of fNIRS Signals into Functional and Systemic Components Based
% on Differences in Hemodynamic Modalities. PLoS One, 7, e50271.
a		= 1/(kf-ks);
b		= [-ks 1; -kf*ks kf];
Hb		= [dHbO; dHbR];
F		= a*b*Hb;
% HbOf	= F(1,:);
% HbRf	= F(2,:);
a = 1/(kf-ks);
b = [kf 1; kf*ks -ks];
S = a*b*Hb;
% dHbOs = S(1,:);
% dHbRs = S(2,:);
