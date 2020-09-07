classdef lsldert_pubclient < lsldert_abstract_client
    properties
        default_port = 5556;
        default_hostname = 'lsldert00.local';
        context;
        socket;
        hostname;
        hostip;
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
            this.hostip=gethostbyname(this.hostname);
            if isempty(this.hostip)
                ME = MException('lsldert_pubclient:hostname_lookup_failure',...
                    'cannot resolve IP address for host %s ', hostname);
                throw(ME);
            end
            % zmq workaround to prevent infinite wait:
            % probe if hostname:port is available
            con=pnet('tcpconnect',this.hostip,port);
            if con < 0
                ME = MException('pupil_remote_control:connect_error',...
                    'cannot connect to tcp://%s:%d, check hostname, port number and/or port number in pupil-capture remote control plugin',...
                    this.hostip, port);
                throw(ME);
            end
            pnet(con,'close');
            
            this.context=ZMQ.context(1);
            this.socket=this.context.socket(ZMQ.PUB);
            uri=sprintf('tcp://%s:%d',this.hostip,port);
            this.socket.connect(uri);
            this.socket.getReceiveTimeOut();

            wait_for_connection(60);
            
            function wait_for_connection(seconds)
                subsuri=sprintf('tcp://%s:%d',this.hostip,port+1);  % assume port+1
                subs=this.context.socket(ZMQ.SUB);
                subs.connect(subsuri);
                timeout=100;
                subs.setReceiveTimeOut(timeout);
                topic='INIT';
                subs.subscribe(topic);
                % try to get a INIT response back before proceeding
                for ii=1:seconds*1000/timeout
                    str=sprintf('%s %d',topic, ii);
                    this.send(str);
                    result=subs.recvStr();
                    result=erase(char(result),char(0));
                    if strcmp(result,str)
                        return
                    end
                end
                fprintf('received ''%s''\n',result);
                error('no response from proxy %s\n', subsuri);
            end
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
