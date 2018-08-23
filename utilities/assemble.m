function [a,u] = assemble(val,x,varargin)
% [A,U] = ASSEMBLE(Y,X)
%
% gathers all Y values for each unique value U in X, and puts the
% mean in A.
%
% [A,U] = ASSEMBLE(Y,X,'fun',@FUNCTION)
%
% Performs function @FUNCTION on the assembled Y values, rather than the
% default mean (e.g. @sum, @std)
%
% Note: this is a really simple function, but this helps me type less...
%
%
% See also ACCUMARRAY, UNIQUE, AVENGERS ASSEMBLE

fun				= keyval('fun',varargin,@mean); % function


[m,n] = size(val);
if n>m
	val = transpose(val);
end
n		= size(val,1);
nx		= size(x,2);
if nx==n
	x = transpose(x);
end


[u,~,subs]		= unique(x,'rows');
a				= accumarray(subs,val,[],fun);
