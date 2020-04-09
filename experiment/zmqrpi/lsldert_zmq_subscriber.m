classdef lsldert_zmq_subscriber < handle
    properties (Access=protected)
        default_port = 5556;
        default_hostname = 'raspi4';
        context;
        socket;
        hostname;
        tmax = Inf;
    end
    
    methods (Access=public)
        function this = lsldert_zmq_subscriber(hostname,port)
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
            
%             % zmq workaround to prevent infinite wait:
%             % probe if hostname:port is available
%             con=pnet('tcpconnect',hostname,port);
%             if con < 0
%                 ME = MException('pupil_remote_control:connect_error',...
%                     'cannot connect to tcp://%s:%d, check hostname, port number and/or port number in pupil-capture remote control plugin',...
%                     hostname, port);
%                 throw(ME);
%             end
%             pnet(con,'close');
%             
            this.context=ZMQ.context(1);
            this.socket=this.context.socket(ZMQ.SUB);
            uri=sprintf('tcp://%s:%d',hostname,port);
            this.socket.connect(uri);
            this.socket.setReceiveTimeOut(1000);
        end
        
        function result=recv(this)
            % SEND send a message string msg to an lsldert client
            % result = obj.SEND(msg [, extra_args...]) sends  message string
            % msg to lsldert. The extra arguments can be text or numeric.
            % all numeric data is sent as a row-major (lexicographically)
            % ordered vector. This is the common
            % ordering for programming languages like C, C++,
            % Python (NumPy) and Pascal.
            import org.zeromq.ZMQ;
            result=this.zmq_recv_multi();
        end
        
        function result=recv_str(this, dontwait)
            import org.zeromq.ZMQ;
            if nargin < 2
                dontwait=true;
            end
            if dontwait
                flags=ZMQ.DONTWAIT;
            else
                flags=0;
            end
            result=this.socket.recvStr(flags);
        end
        
        function result=zmq_recv_multi(this)
            import org.zeromq.ZMQ;
    %TODO read multiple frames
            result=this.socket.recv();
        end
    end
    
end