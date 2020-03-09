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
        end
        
        function result=send(this, msg)
            % SEND send a message string to an lsldert server
            % result = obj.SEND(msg) sends  message string str to lsldert
            % and returns the  message string returned from the remote
            % side.
            tstart=tic;
            this.socket.send(sprintf('%s\0',msg));
            
            %rbytes=this.socket.recv();
            telapsed=toc(tstart);
            if telapsed > this.tmax
                warning('lsldert_client.send: round trip time exceeded max, %2.1f>%2.1f ms', 1e3*telapsed, 1e3*this.tmax);
            end
            % telapsed
            result='';
        end
        
    end
end
