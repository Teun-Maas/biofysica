classdef lsldert_pubclient < lsldert_abstract_client
    properties
        default_port = 5556;
        default_hostname = 'lsldert00.local';
        context;
        socket;
        hostname;
        tmax = Inf;
        
    end
    
    methods
        function this = lsldert_pubclient(hostname,port)
            % LSLDERT_CLIENT class constructor
            %
            % obj = lsldert_client('lslsdert-host.local');
            % creates a connection to the lsldert server listening
            % on port 5555 at remote host 'lslsdert-host.local'.
            % the port number can optionally be specified as a 2nd argument
            % obj = lsldert_client('lsldert-host.local', 5555);
            
            import org.zeromq.ZMQ;
            if nargin < 2
                port = this.default_port;
            end
            if nargin < 1
                hostname = this.default_hostname;
            end
            this.hostname=hostname;
            
            % zmq workaround to prevent infinite wait:
            % probe if hostname:port is available
            con=pnet('tcpconnect',hostname,port);
            if con < 0
                ME = MException('pupil_remote_control:connect_error',...
                    'cannot connect to tcp://%s:%d, check hostname, port number and/or port number in pupil-capture remote control plugin',...
                    hostname, port);
                throw(ME);
            end
            pnet(con,'close');
            
            this.context=ZMQ.context(1);
            this.socket=this.context.socket(ZMQ.PUB);
            uri=sprintf('tcp://%s:%d',hostname,port);
            this.socket.connect(uri);
            this.socket.getReceiveTimeOut();

%           function wait_for_connection
%              subsuri=sprintf('tcp://%s:%d',hostname,port+1);  % assume port+1
%              subs=this.context.socket(ZMQ.SUB);
%              subs.connect(subsuri);
%              subs.setReceiveTimeout(1000);
%              topic='INIT';
%              subs.subscribe(topic);
%              for ii=1:10
%                 str=sprintf('%s %d',topic, ii);
%                 this.send(str);
%                 result=subs.recv;
%                 fprintf('received %s\n',result);
%                 if strcmp(result,str)
%                    return
%                 end
%              end
%           end
        end
        
        function [result,telapsed]=send(this, msg, varargin)
            % SEND send a message string msg to an lsldert client
            % result = obj.SEND(msg [, extra_args...]) sends  message string
            % msg to lsldert. The extra arguments can be text or numeric.
            % all numeric data is sent as a row-major (lexicographically)
            % ordered vector. This is the common
            % ordering for programming languages like C, C++,
            % Python (NumPy) and Pascal.
            tstart=tic;
            this.zmq_write_multi(this.socket, msg, varargin{:});
            telapsed=toc(tstart);
            if telapsed > this.tmax
                warning('lsldert_client.send: round trip time exceeded max, %2.1f>%2.1f ms', 1e3*telapsed, 1e3*this.tmax);
            end
            result='';
        end
        
    end
end
