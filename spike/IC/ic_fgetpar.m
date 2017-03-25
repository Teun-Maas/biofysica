function param = ic_fgetpar(fp, par, format)

% FUNCTION param = fgetpar(fp, par, format)
%
%   read parameter(s) "par" form text file fp   
%   the search procedure for parameter names
%   is case sensitive.
%   
%   Jeroen Goossens


% search from beginning 
frewind(fp);
par = ['#' par];

while ~feof(fp)
  % search par identifier 
  ident=fscanf(fp,'%s',1);
  if strcmp(ident,par)==1
    % find '=' sign 
    symb=fscanf(fp,'%s',1);
    if strcmp(symb,'=')==1
      % read value(s)
      param = fscanf(fp,format,1);
      return
    end
    return
  end
end

if ~exist('param', 'var'), param=[]; end
