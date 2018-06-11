function l=lsl_list(varargin)
% LSL_LIST - resolves lsl streams on the local network and returns them in
% a Nx2 cell array.
% To resolve streams one of three calling methods can be used:
% lsl_list() uses lsl_resolve_all(lib);
% lsl_list(predicate_string) uses lsl_resolve_bypred(lib,predicate_string);
% lsl_list(property,value) uses lsl_resolve_byprop(lib,property,value);
%
% Examples:
% lsl_list(); % resolve all streams
% lsl_list('name','Digital Events 0');
% lsl_list('type','Digital Events @ clockpi');
% lsl_list('type=''Digital Events @ clockpi'' and name=''Digital Events 0''')
%
% For more information on predicate strings see the help on lsl_resolve_by_pred.
%
% SEE ALSO: LSL_RESOLVE_ALL, LSL_RESOLVE_BYPRED, LSL_RESOLVE_BYPROP
  info=lsl_resolver(varargin{:});
  l=info.list();
  delete(info);
end
