function h = pa_errorbar3(bardata,errdata)
% H = PA_ERRORBAR3(BARDATA,ERRDATA)
%
% draw a bar3 graph of data BARDATA, with errorbars ERRDATA
%
% see also BAR3, ERRORBAR

% (c) 2011-04-27 Marc van Wanrooij

[Ngroups,Nbars] = size(bardata)  %#ok<ASGLU>
h				= bar3(bardata,1);
hold on;
% h = get(h(Ngroups),'children')

for ii = 1:Ngroups
	x = get(h(ii),'XData');
	y = get(h(ii),'YData');
	z = get(h(ii),'ZData');
	
% 	whos x y z
	x = mean([min(x(:)) max(x(:))]);
	
	n = round(length(y)/Nbars);
	
	plot3()
% 	x = mean(x(:,[1 3]))
% 	y = mean(y(:,[1 3]))
	
% 	if ~isempty(l)
% 		for jj = 1:Nbars
% 			[ii jj]
% 					l(jj)
% 
% % 			k = get(l(jj))
% %% 			if ~isempty(k)
% % 				x = get(k)
% % 			end
% 		end
% 	end
% 	x = get(get(h(ii),'children'),'xdata');
% 	x = mean(x([1 3],:));
% 	y = bardata(:,ii);
% 	e = errdata(:,ii);
% 	errorbar(x, y, e, 'k', 'linestyle', 'none');
end
