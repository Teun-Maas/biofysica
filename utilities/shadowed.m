function shipped = shadowed(this)
% SHIPPED = SHADOWED(THIS)
%
% stores a handle to the shipped MATLAB funtion of the same name%
%
% https://nl.mathworks.com/matlabcentral/answers/149059-how-do-i-invoke-a-shadowed-core-matlab-function-not-built-in-from-an-overloaded-function-of-same

this = fcheckext(this,'.m');
% persistent shipped
% if isempty(shipped)
list = which(this, '-all'); % find all the functions which shadow it
f = strncmp(list, matlabroot, length(matlabroot)); % locate 1st in list under matlabroot
list = list{find(f, 1)}; % extract from list the exact function we want to be able to call
here = cd(list(1:end-length(this))); % temporarily switch to the containing folder
shipped = str2func(this(1:end-2)); % grab a handle to the function
cd(here); % go back to where we came from
% end