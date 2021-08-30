classdef lsl_resolver < handle
    % LSL_RESOLVER - the lsl_resolver class uses the lsl_resolve_* functions to
    % find (lists) of clients, that can easily be referenced by indexing.
    %
    %Examples:
    %info=lsl_resolver('type=''Digital Events @ clockpi'' and name=''Digital Events 0''');
    %info=lsl_resolver('type=''Pupil Capture @ dcn-eyebrain'' and name=''Pupil Primitive Data - Eye 0''')
    %info=lsl_resolver('type=''Pupil Gaze @ LAPTOP-CBLJOAUJ''');
    %
    %info is returned as a cell array of LSL_STREAMINFO objects that can be used to creat
    % LSL_ISTREAM and LSL_OSTREAM obects.
    %
    % See also: LIST, LIST_VERBOSE, LSL_SESSION, LSL_ISTREAM, LSL_OSTREAM
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
            % LIST - List basic information of all streams found in an array of structures
            % with name and type of the streams
            %
            %   l=info.list();
            %   if isempty(l)
            %       error('no streams found');
            %   end
            %
            %   for i=1:size(l,1)
            %       fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
            %   end
            %
            %   n=input('enter stream number to acquire: ');
            %   the_stream=lsl_istream(info{n});
            
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
            % LIST_VERBOSE - same as lsl_resolver.list, but with text output
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
