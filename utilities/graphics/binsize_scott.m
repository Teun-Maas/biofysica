function [N,xi,h,k] = binsize_scott(x,varargin)
% [H
%
% See also: https://en.wikipedia.org/wiki/Histogram


%% Initialization
disp = keyval('disp',varargin,false);

%% Bin width / number of bins
sigma	= std(x); % sample standard deviation
n		= numel(x); % number of elements
h		= 3.5*sigma/(n^(1/3)); % bin width
k		= ceil((max(x)-min(x))/h); % number of bins

%% Histogram
xi		= linspace(min(x),max(x),k);
N		= hist(x,xi);

%% Graphics
if disp
stairs(xi,N)
axis square;

end