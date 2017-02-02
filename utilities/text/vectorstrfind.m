function v = vectorstrfind(C,str,varargin)
% V = VECTORSTRFIND(TEXTCELLARRAY,PATTERN)
%
% Returns the cell indices of any occurrences of the string PATTERN in
% TEXTCELLARRAY. 
%
% Optional:
% V = VECTORSTRFIND(TEXTCELLARRAY,PATTERN,'logic',true)
%
% returns logical vector
%
% See also STRFIND

logicFlag = keyval('logic',varargin,false);

idx		= strfind(C,str);

v		= not(cellfun('isempty', idx));

if ~logicFlag
	v = find(v);
end