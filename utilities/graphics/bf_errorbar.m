function h = pa_errorbar(bardata,errdata,varargin)
% H = PA_ERRORBAR(BARDATA,ERRDATA)
%
% draw a bar graph of data BARDATA, with errorbars ERRDATA
%
% see also BAR, ERRORBAR

% 2013-01-04 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Optional, undocumented:
col				= pa_keyval('Color',varargin);
xlabels			= pa_keyval('xlabel',varargin);
bsize			= pa_keyval('bsize',varargin);

[Ngroups,Nbars] = size(bardata); %#ok<ASGLU> % Number of data points 
[~,Nerr]		= size(errdata); % determine: one error, or upper and lower errors
h				= bar(bardata,bsize);
if ~isempty(col)
	set(h,'FaceColor',col);
end
hold on;
for ii = 1:Nbars
	x = get(get(h(ii),'children'),'xdata');
	x = mean(x([1 3],:));
	y = bardata(:,ii);
	indx = (1:2)+(ii-1)*2;
	if (Nerr/Nbars)==2
		e = errdata(:,indx);
		errorbar(x, y, e(:,1),e(:,2), 'k', 'linestyle', 'none');
	else
		e = errdata(:,ii);
		errorbar(x, y, e, 'k', 'linestyle', 'none');
	end
end

%% X-labels
if ~isempty(xlabels)
	n = numel(xlabels);
	set(gca,'XTick',1:n,'XTickLabel',xlabels);
end