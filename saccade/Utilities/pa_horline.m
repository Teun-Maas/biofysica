function h = pa_horline(y, style)
% PA_HORLINE(Y)
%
% Plot horizontal line through current axis at height Y.
%
% PA_HORLINE(...,'LineSpec') uses the color and linestyle specified by 
% the string 'LineSpec'. See PLOT for possibilities.
%
% H = PA_HORLINE(...) returns a vector of lineseries handles in H.
%
% See also PA_VERLINE, PA_UNITYLINE
%

% (c) 2011 Marc van Wanrooij

%% Initialization
if nargin < 2, style = 'k--'; end
if nargin < 1, y = 0; end

y       = y(:)'; % Create a column vector
n_as    = get(gca,'Nextplot');
n_fi    = get(gcf,'Nextplot');
oldhold = ishold;
hold on
x_lim   = get(gca,'XLim');
for i = 1:length(y)
    h       = plot(x_lim, [y(i);y(i)], style);
end
set(gca,'Nextplot',n_as);
set(gcf,'Nextplot',n_fi);
if ~oldhold
	hold off;
end

