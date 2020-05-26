classdef lsl_resolver < handle
    properties
        infos
    end
    
    methods
        function this=lsl_resolver(arg1, arg2)
            lib=lsl_loadlib();
            if nargin==0
                this.infos=lsl_resolve_all(lib);
            elseif nargin==1
                this.infos=lsl_resolve_bypred(lib,arg1);
            else
                this.infos=lsl_resolve_byprop(lib,arg1,arg2);
            end
        end

        function l=list(this)
            n=length(this.infos);
            if n==0
                l=[];
                return
            end
            l=struct();
            for i=1:n
                l(i).name=this.infos{i}.name;
                l(i).type=this.infos{i}.type;
            end
        end
        
        function l=list_verbose(this)
            l=this.list();
            for i=1:numel(l)
                disp(l(i));
            end
        end
        
        function varargout=subsref(this,s)
            switch s(1).type
                case {'.','()'}
                    varargout = {builtin('subsref',this,s)};
                    
                case '{}'
                    if length(s) == 1
                        % Implement obj{indices}
                        varargout={this.infos{s(1).subs{1}}};
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName
%                         ind=s(1).subs{1}
%                         nam=s(2).subs
                        varargout={this.infos{s(1).subs{1}}.(s(2).subs)};
                    else
                        % Use built-in for any other expression
                        varargout = {builtin('subsref',this,s)};
                    end
                otherwise
                    error('Not a valid indexing expression')
                    
            end
            
        end
    end
end
