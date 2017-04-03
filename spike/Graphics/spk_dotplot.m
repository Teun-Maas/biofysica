function h = pa_spk_dotplot(spikes,varargin)
% PA_SPK_DOTPLOT(SPK)
%   Plot a spike train raster of the spikes in SPK. 
% SPK is a MxN matrix 
%, where each row is a row
%   in the dot display.
%
% See also PA_SPK_RASTERPLOT

% (c) 2011 Marc van Wanrooij
% e-mail: m.vanwanrooij@gmail.com

marker = pa_keyval('marker',varargin);
if isempty(marker)
	marker = '.';
end
markersize = pa_keyval('markersize',varargin);
if isempty(markersize)
	markersize = 1;
end

hold on;
N	= size(spikes,1);
for ii=1:N,
   Trace	= find(spikes(ii,:));			% Find the indices of the samples that contain a spike
   Row		= repmat(ii,size(Trace));		% Create a vector the same size as Trace, with Trial Number as only value in this vector
   plot(Trace,Row,['k' marker], 'markersize',markersize);	% Now plot at each time and trial number where a spike occurs a single dot
end;
% Labeling of the graph
box on;
xlabel('Time (ms)');
ylabel('Trial #');
ylim([0 N+1]);

% Check for output
if nargout, 
	h = gca; % If output is asked, return the current axis handle
end; 
