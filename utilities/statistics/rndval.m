function val = rndval(a,b,n)
% RNDVAL Random uniform values from the interval [A,B]
%
% VAL = RNDVAL(A,B,N)
% Obtain N uniform values from the interval [A,B]. 
%
% Example
% To create a random distribution for fixation LED offsets during an
% experiment to avoid expectations from subjects, choose an interval of
% 300-800 msec:
%   ledoff = rndval(300,800,[50 2]);
%
% See also RAND

% Marc van Wanrooij, 2007

b	= b+1;
val = floor(a+(b-a)*rand(n));

