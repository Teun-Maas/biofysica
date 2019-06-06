function h = verline(x, style)
% VERLINE(Y)
%
% Plot vertical line through current axis at width X.
%
% VERLINE(...,'LineSpec') uses the color and linestyle specified by 
% the string 'LineSpec'. See PLOT for possibilities.
%
% H = VERLINE(...) returns a vector of lineseries handles in H.
%
% See also HORLINE, UNITYLINE
%

% (c) 2011 Marc van Wanrooij

%% Initialization
if nargin < 2, style = 'k--'; end
if nargin < 1, x = 0; end


x       = x(:)'; % Create a column vector
n_as    = get(gca,'Nextplot');
n_fi    = get(gcf,'Nextplot');
oldhold = ishold;
hold on;
y_lim   = get(gca,'YLim');
for i = 1:length(x)
    h       = plot([x(i);x(i)], y_lim, style);
end
set(gca,'Nextplot',n_as);
set(gcf,'Nextplot',n_fi);
if ~oldhold
	hold off;
end

if nargout==0
	clear h
end