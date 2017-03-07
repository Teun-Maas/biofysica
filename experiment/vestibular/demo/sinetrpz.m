function [ y ] = sinetrpz(A, T, t, Ts1, Ts2)
%MERFELD Produce a sine function tapered by a trapezoid
%   merfeld2(A, T, t, Ts1, Ts2) produces a sine function
%   y=A*sin(2*pi*t./T) that is windowed using a trapezoid
%   of base 0 and height 1 on the interval [0 Ts1 max(t)-Ts2 max(t)]
%   Hence the starting rise time  is Ts1, the stopping rise time is Ts2.

if nargin  < 4
    Ts1 = 0;
end
if nargin < 5
    Ts2 = 0;
end

y=A*sin(2*pi*t./T);
w=l_trapmf(t,[0 Ts1 max(t)-Ts2 max(t)]);
y=y.*w;

end


% local 'trapmf' function
function y = l_trapmf(x, params)
%l_trapmf Trapezoidal membership function.
%   l_trapmf(X, PARAMS) returns a matrix which is the trapezoidal
%   membership function evaluated at X. PARAMS = [A B C D] is a 4-element
%   vector that determines the break points of this membership function.
%   We require that A <= B and C <= D. If B >= C, this membership
%   function becomes a triangular membership function that could have
%   a height less than unity. (See the example below.)
%
%   For example:
%
%       x = (0:0.1:10)';
%       y1 = l_trapmf(x, [2 3 7 9]);
%       y2 = l_trapmf(x, [3 4 6 8]);
%       y3 = l_trapmf(x, [4 5 5 7]);
%       y4 = l_trapmf(x, [5 6 4 6]);
%       plot(x, [y1 y2 y3 y4]);
%       set(gcf, 'name', 'l_trapmf', 'numbertitle', 'off');
%
%   See also DSIGMF, EVALMF, GAUSS2MF, GAUSSMF, GBELLMF, MF2MF, PIMF, PSIGMF,
%   SIGMF, SMF, TRIMF, ZMF.

%   Roger Jang, 6-28-93, 10-5-93, 4-14-94.
%   Copyright 1994-2002 The MathWorks, Inc.
%   $Revision: 1.22 $  $Date: 2002/04/14 22:21:13 $
%   G?nter Windau, 20170306, syntax corrections


if nargin ~= 2
    error('Two arguments are required by the trapezoidal MF.');
elseif length(params) < 4
    error('The trapezoidal MF needs at least four parameters.');
end

a = params(1); b = params(2); c = params(3); d = params(4);

if a > b
    error('Illegal parameter condition: a > b');
elseif c > d
    error('Illegal parameter condition: c > d');
end

y1 = zeros(size(x));
y2 = zeros(size(x));

% Compute y1
index = find(x >= b);
if ~isempty(index)
    y1(index) = ones(size(index));
end
index = find(x < a);
if ~isempty(index)
    y1(index) = zeros(size(index));
end
index = find(a <= x & x < b);
if ~isempty(index) && a ~= b
    y1(index) = (x(index)-a)/(b-a);
end

% Compute y2
index = find(x <= c);
if ~isempty(index)
    y2(index) = ones(size(index));
end
index = find(x > d);
if ~isempty(index)
    y2(index) = zeros(size(index));
end
index = find(c < x & x <= d);
if ~isempty(index) && c ~= d
    y2(index) = (d-x(index))/(d-c);
end

% Compute y
y = min(y1, y2);

end