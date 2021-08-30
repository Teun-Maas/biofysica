function y = startsWith(string, pattern)
% INPUTS:
%   - string : char-array
%   - pattern : prefix-pattern to find in string
%
% OUTPUTS
%   - y : true is pattern was found in the beginning 

y = false;
o = strfind(string, pattern);
if length(o) == 1 && o == 1
   y = true;
end