function plotfart(mc)
% PLOTFART plot the double polar-limits of the FART-setup
%
% PLOTFART(MC)
% Plot the limits with color MC.
%
% Example
%   plotfart('r');
%
% See also GENEXPERIMENT, PLOTGRID
%

% Author  Marc van Wanrooij
% Copyright 2007

if nargin<1
    mc = [0.5 0.5 0.5];
end
oldhold = ishold;
if ~oldhold
    hold on;
end
axis([-100 100 -100 100]);
% for i = 1:(29+29)
%     h=horline(((i-1)*2.5)-57.5,'k-');set(h,'Color',[0.7 0.7 0.7]);
% end

h = plot([-90 -5],[0 85],'k-',...
    [5 90],[85 0],'k-',...
    [32.5 90],[-57.5 0],'k-',...
    [-90 -32.5],[0 -57.5],'k-',...
    [-32.5 32.5],[-57.5 -57.5],'k-',...
    [-5 5],[85 85],'k-'); set(h,'LineWidth',2,'Color',mc);
fixdist = 8;
elfix       = -57.5:0.1:85;
azfix       = asind(sind(fixdist).*cosd(elfix));
h           = plot(azfix,elfix,'k-');set(h,'LineWidth',2,'Color',mc);
elfix       = -57.5:0.1:85;
azfix       = asind(sind(-fixdist).*cosd(elfix));
h           = plot(azfix,elfix,'k-');set(h,'LineWidth',2,'Color',mc);

axis square;
box on;
if ~oldhold
    hold off;
end