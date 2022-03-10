classdef lslder_remote_control < handle
    % LSLDER_REMOTE_CONTROL is a class to communicate with a zmqrpi remote trigger
    % program
    % See also: ZMQ_RPI_REMOTE_CONTROL, PUPIL_REMOTE_CONTROL, SEND
    %
    % LSLDER_REMOTE_CONTROL uses JeroMQ, a pure java implementation of
    % libzmq. https://github.com/zeromq/jeromq

    % Version 1.1 GW/20220310-1
    %
    
    properties (Access=protected)
        default_port = 50020;
        context;
        socket;
        hostname;
        tmax = Inf;
    end
    
    methods
        function this = lslder_remote_control(hostname, port)
            % PUPIL_REMOTE_CONTROL class constructor
            % 
            % obj = lslder_remote_control('pupil-host.local');
            % creates a connection to the Pupil Remote plugin listening
            % on port 50020 at remote host 'pupil-host.local'.
            % the port number can optionally be specified as a 2nd argument
            % obj = lslder_remote_control('pupil-host', 50021);
            
            import org.zeromq.ZMQ;
            if nargin < 2
                port = this.default_port;
            end
            this.hostname=hostname;
            
            % zmq workaround to prevent infinite wait:
            % probe if hostname:port is available 
            con=pnet('tcpconnect',hostname,port);
            if con < 0
                ME = MException('lslder_remote_control:connect_error',...
                    'cannot connect to tcp://%s:%d, check hostname, port number and/or port number in pupil-capture remote control plugin',...
                    hostname, port);
                throw(ME);
            end
            pnet(con,'close');
            
            this.context=ZMQ.context(1);
            this.socket=this.context.socket(ZMQ.REQ);
            uri=sprintf('tcp://%s:%d',hostname,port);
            this.socket.connect(uri);
            this.socket.getReceiveTimeOut();
        end
        
        function delete(this)
            this.socket.close();
            this.context.term();
        end
        
        function result=send(this, msg)
            % SEND send a message string to Pupil Remote
            % result = obj.SEND(msg) sends  message string str to Pupil Remote 
            % and returns the  message string returned from the remote
            % side.
            tstart=tic;
            this.socket.send(sprintf('%s\0',msg));

            rbytes=this.socket.recv();
            telapsed=toc(tstart);
            if telapsed > this.tmax
                warning('lslder_remote_control: round trip time exceeded max, %2.1f>%2.1f ms', 1e3*telapsed, 1e3*this.tmax);
            end
            % telapsed
            result=reshape(char(rbytes),1,[]);
        end
        
        function result = pulseIR(this, mask, duration, marker)
            if nargin < 2
                mask = 3;
            end
            if nargin < 3
                duration = 0.5;
            end
            if nargin < 4
                marker = sprintf("PULSEIR_%d", mask);
            end
            str = sprintf("I %d %d %s", mask, duration*1000, marker);
            result = this.send(str);
        end
        
        function result = beep(this, freq, duration, marker)
            % BEEP - play a beep on the raspberry
            % result = obj.beep(freq, duration);  
            % result = obj.beep(880, 0.6); % plays a 880 Hz tone for 0.6
            % seconds
            % If freq or duration are omitted, default values are 500Hz and
            % 1 second.
            if nargin < 2
                freq = 440;
            end
            if nargin < 3
                duration = 0.5;
            end
            if nargin < 4
                marker = sprintf("BEEP_%d",freq);
            end
            str = sprintf("B %d %d %s", freq, duration*1000, marker);
            result = this.send(str);
        end
        
        function result = digitalout(this, value, marker)
            % DIGITALOUT - set digital output to high or low level
            if nargin < 3
                marker = sprintf("DOUT_%d", value);
            end
            str = sprintf("D %d %s", value, marker);
            result = this.send(str);
        end
        
        function set_digitalin_marker(this, marker)
            % SET_DIGITALIN_MARKER - set the marker name that is sent on digital in events
            str = sprintf("M %s", marker);
            result = this.send(str);
        end
    end
end
