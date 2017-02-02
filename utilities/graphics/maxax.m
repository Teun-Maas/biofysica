function ax = pa_maxax
%            AX = PA_MAXAX
%
% Set axis to maximal absolute values.

% (c) 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

ax  = axis;
mx  = max(abs(ax));
axis([-mx mx -mx mx]);
ax = axis;
