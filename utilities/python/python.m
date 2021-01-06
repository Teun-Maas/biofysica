function [result, status] = python(varargin)
%   python Execute Python command and return the result.
%   python(PYTHONFILE) calls python script specified by the file PYTHONFILE
%   using appropriate python executable.
%
%   python(PYTHONFILE,ARG1,ARG2,...) passes the arguments ARG1,ARG2,...
%   to the python script file PYTHONFILE, and calls it by using appropriate
%   python executable.
%
%   RESULT=python(...) outputs the result of attempted python call.  If the
%   exit status of python is not zero, an error will be returned.
%
%   [RESULT,STATUS] = python(...) outputs the result of the python call, and
%   also saves its exit status into variable STATUS.
%
%   If the Python executable is not available, it can be downloaded from:
%     http://www.python.org
%
%   See also system, java, mex, python
global BIOFYSICA_PY_EXE
if isempty(BIOFYSICA_PY_EXE)
    error('BIOFYSICA:python:NoExecutable. Check that biofysica toolbox has correctly defined the python interpreter');
end

if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

cmdString = '';

% Add input to arguments to operating system command to be executed.
% (If an argument refers to a file on the MATLAB path, use full file path.)
for i = 1:nargin
    thisArg = varargin{i};
    if ~ischar(thisArg)
        error('BIOFYSICA:python:InputsMustBeStrings');
    end
    if i==1
        if exist(thisArg, 'file')==2
            % This is a valid file on the MATLAB path
            if isempty(dir(thisArg))
                % Not complete file specification
                % - file is not in current directory
                % - OR filename specified without extension
                % ==> get full file path
                thisArg = which(thisArg);
            end
        else
            % First input argument is PythonFile - it must be a valid file
            error('BIOFYSICA:python:FileNotFound %s', thisArg);
        end
    end
    
    % Wrap thisArg in double quotes if it contains spaces
    if isempty(thisArg) || any(thisArg == ' ')
        thisArg = ['"', thisArg, '"']; %#ok<AGROW>
    end
    
    % Add argument to command string
    cmdString = [cmdString, ' ', thisArg]; %#ok<AGROW>
end

% Check that the command string is not empty
if isempty(cmdString)
    error('BIOFYSICA:python:NoPythonCommand %s',cmdString);
end

cmdString = [BIOFYSICA_PY_EXE, ' ', cmdString];
cmdString
[status, result] = system(cmdString);

if nargout < 2 && status~=0
    error('MATLAB:python:ExecutionError, %s, %s', result, cmdString);
end

end
