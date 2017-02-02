% function nirs_detrend(data,varargin)
function data = nirs_detrend(data,varargin)

% DATA = NIRS_RCS(DATA)
%
% Perform polynomial detrend subtraction on DATA
%
% See also NIRS_OPTODELAYOUT, NIRS_OPTODEDISTANCE

%% Initialization

if nargin<1
	close all
	gr = 'CI';
	d				= keyval('dir',varargin,['/Users/marcw/DATA/Roos Cartignij/NIRS sessie' filesep gr]);
	cd(d);
	
	%% Do stuff
	expdir			= dir;
	expnames		= {expdir.name};
	idxO2				= vectorstrfind(expnames,'.');
	expnames(idxO2)	= [];
	expname = expnames{5};
	cd(expname)
	fname	= dir('rcs*.mat');
	fname	= fname.name;
	load(fname)
	data = data_rcs;
end;

dispFlag	= keyval('disp',varargin,false);
method		= keyval('method',varargin,'closest');

trial = data.trial{1};
x = data.time{1};
ntrial = size(trial,1);
for ii = 1:ntrial
	y = trial(ii,:);
	sel = ~isnan(x);
	p			= polyfit(x(sel),y(sel),20);
	a			= y-polyval(p,x);
	trial(ii,:) = a;
end
data.trial{1} = trial;