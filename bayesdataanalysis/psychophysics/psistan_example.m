function psistan_example

%% Simulation
[x,y,s] = getdata;

% keyboard
close all;
cd /Users/marcw/Documents/MATLAB


		centroid = 'median';
		guessrate = 0;
		lapserate = 0;
		fun = @logisticfun;

t = tic;
samples = psistan(x,y,s,...
	'gamma',guessrate,...
	'lambda',lapserate,...
	'function',fun,...
	'numSavedSteps',5000,...
	'showCentroid',centroid,...
	'showDiag',false);
e = toc(t)
%%


%%
% savegraph(mfilename,'png');
keyboard

%%
function [x,y,s] = getdata

%% Bayesian Causal Modeling workshop Michael Lee & Eric
% load data
p = which('data_x.txt') % does it exist?
[pathstr,~,~] = fileparts(p);
fname	= 'data_x.txt';
x		= dlmread(fullfile(pathstr,fname),'\t');
fname	= 'data_n.txt';
n		= dlmread(fullfile(pathstr,fname),'\t');
fname	= 'data_r.txt';
r		= dlmread(fullfile(pathstr,fname),'\t');
nsubjs	= size(x,1);
% 'transform' data
S		= repmat(transpose(1:nsubjs),1,size(x,2)); % subject vector
X		= x(:);
Y		= [r(:) n(:)]; % 2-dimensional y containing rate and number of trials
S		= S(:);
% remove NaNs
sel		= ~isnan(x);
x		= X(sel);
y		= Y(sel,:);
s		= S(sel,:);

