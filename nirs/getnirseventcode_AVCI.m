function [t,utxt,ustim,stimV,stimA,stimAV] = getnirseventcode_AVCI(dname,fname)

if nargin<1
	dname = '/Users/marcw/DATA/Roos Cartignij/NIRS sessie/CI/CI_AV_6';
end
cd(dname)

if nargin<2
	f			= dir('*order.xlsx');
	fname	= f.name;
end
N = xlsread(fname);

stim = N(:,1);
% xx
% Convention
% 0, 1, 2, 3 = 0%, 25%, 50%, 75%
% decades = Visual, units = auditory

stimV			= floor(stim/10);
stimA			= mod(stim,10);
stimAV			= stimA & stimV;
stim			= [stimV stimA];
[ustim,~,t]		= unique(stim,'rows');
nstim			= size(ustim,1);

%% Text
prctxt = {'00','25','50','75'};

for ii = 1:nstim
	utxt(ii,:) = ['V' prctxt{ustim(ii,1)+1} 'A' prctxt{ustim(ii,2)+1}]; %#ok<AGROW>
end


% t is from 1 to 11
% 1 = V0, A25
% 2 = V0, A50
% 3 = V0, A75
% 4 = V25, A0
% 5 = V25, A25
% 6 = V25, A75
% 7 = V50, A0
% 8 = V50, A50
% 9 = V75, A0
% 10 = V75, A25
% 11 = V75, A75
