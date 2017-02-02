function h = pa_revline(style)
% PA_REVLINE
%
% Plot diagonal reverse unity line through current axis.
%
% PA_REVLINE('LineSpec') uses the color and linestyle specified by 
% the string 'LineSpec'. See PLOT for possibilities.
%
% See also PA_HORLINE, PA_VERLINE, PA_UNITYLINE
%

% (c) 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com
%% Initialization
if nargin < 1, style='k--'; end
ho = ishold;
if ~ho
    hold on
end

%% Set axis
ax		= axis;
if ax(1)~=ax(3)
    min_x = min(ax([1 3]));
    ax([1 3]) = min_x;
end
if ax(2)~=ax(4)
    max_x = max(ax([2 4]));
    ax([2 4]) = max_x;
end
% axis(ax);
ax([1 2]) = ax([2 1]);

%% Plot unityline
h = plot([ax(1) ax(2)],[ax(3) ax(4)],style);
% axis(ax); % reset axis for if PLOT-command changed axis

%% Checkout
if ~ho
    hold off
end
    
% 
% if nargin < 1, style='k--'; end
% ho = ishold;
% if ~ho
%     hold on
% end
% ax      = axis;
% % ax
% if ax(1)~=ax(3)
%     min_x       = max(ax([1 3]));
%     ax([1 3])   = min_x;
% end
% if ax(2)~=ax(4)
%     max_x       = min(ax([2 4]));
%     ax([2 4])   = max_x;
% end
% ax([1 2]) = ax([2 1]);
% % axis(ax);
% h=plot([ax(1) ax(2)],[ax(3) ax(4)],style);
% % axis(ax);
% if ~ho
%     hold off
% end
%     