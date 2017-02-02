function h = pa_errorbar2(bardata,errdata,varargin)
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

[Nbars,Ngroups] = size(bardata); %#ok<ASGLU> % Number of data points 
[~,Nerr]		= size(errdata); % determine: one error, or upper and lower errors
h				= bar(bardata,1);
% if ~isempty(col)
% 	set(h,'FaceColor',col);
% end
hold on;
for ii = 1:Ngroups
	ch		= get(h(ii),'children'); % get children of bar group
	if ~isempty(col)
		set(ch,'FaceVertexCData',col{ii});
	end
	x		= get(ch,'xdata'); % get the x position of every bar vertex in the group
	x		= mean(x([1 3],:)); % get the center x  position
	y		= bardata(:,ii); % the y-data
	
	%% plot errorbars
	indx	= (1:2)+(ii-1)*2;
	if (Nerr/Ngroups)==2 % if lower and upper errors
		e	= errdata(:,indx);
		errorbar(x, y, e(:,1),e(:,2), 'k', 'linestyle', 'none');
	else % if one error
		e	= errdata(:,ii);
		errorbar(x, y, e, 'k', 'linestyle', 'none');
	end
end

%% X-labels
if ~isempty(xlabels)
	n = numel(xlabels);
	set(gca,'XTick',1:n,'XTickLabel',xlabels);
end
