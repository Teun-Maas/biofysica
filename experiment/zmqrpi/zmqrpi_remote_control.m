classdef zmqrpi_remote_control < handle
    % ZMQ_RPI_REMOTE_CONTROL is a class to communicate with a zmqrpi remote trigger
    % program
    % See also: PUPIL_REMOTE_CONTROL, SEND
    %
    % ZMQ_RPI_REMOTE_CONTROL uses JeroMQ, a pure java implementation of
    % libzmq. https://github.com/zeromq/jeromq

    % Version 1.0 GW/20191211-1
    %
    
    properties (Access=protected)
        default_port = 50020;
        context;
        socket;
        hostname;
        tmax = Inf;
    end
    
    methods
        function this = zmqrpi_remote_control(hostname, port)
            % PUPIL_REMOTE_CONTROL class constructor
            % 
            % obj = pupil_remote_control('pupil-host.local');
            % creates a connection to the Pupil Remote plugin listening
            % on port 50020 at remote host 'pupil-host.local'.
            % the port number can optionally be specified as a 2nd argument
            % obj = pupil_remote_control('pupil-host', 50021);
            
            import org.zeromq.ZMQ;
            if nargin < 2
                port = this.default_port;
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
            this.socket=this.context.socket(ZMQ.REQ);
            uri=sprintf('tcp://%s:%d',hostname,port);
            this.socket.connect(uri);
            rt=this.socket.getReceiveTimeOut();
rt
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
            this.socket.send(msg);
            rbytes=this.socket.recv();
            telapsed=toc(tstart);
            if telapsed > this.tmax
                warning('pupil_remote_control: round trip time exceeded max, %2.1f>%2.1f ms', 1e3*telapsed, 1e3*this.tmax);
            end
            % telapsed
            result=reshape(char(rbytes),1,[]);
        end
        
        function result = beep(this)
            % BEEP - play a beep on the raspberry
            result = this.send('b');
        end
        
%%% PupilLabs leftovers
        function result = start_recording(this,session_name)
            % START_RECORDING starts a recording session on remote
            % result = obj.start_recording;
            % result = obj.start_recording('my session name');
            if nargin > 1
                msg = sprintf('R %s', session_name);
            else
                msg = 'R';
            end
            result = this.send(msg);
        end
        
        function result = stop_recording(this)
            % STOP_RECORDING stops the current recording session
            % result = obj.stop_recording;
            result = this.send('r');
        end
        
        function result = start_calibration(this)
            % START_CALIBRATION starts the currently selected calibration
            % result = obj.start_calibration;
            result = this.send('C');
        end
        
        function result = stop_calibration(this)
            % STOP_CALIBRATION stops the currently selected calibration
            % result = obj.stop_calibration;
            result = this.send('c');
        end
        
        function result = time_sync(this, reference_time)
            % TIME_SYNC set the remote time stamp clock to reference time
            % result = obj.time_sync(1.234); make time stamps count from
            % 1.234 from now on.
            msg = sprintf('T %18.15f', reference_time);
            result = this.send(msg);
        end
        
        function ts = get_time_stamp(this)
            % GET_TIME_STAMP get remote time stamp
            % ts = obj.get_time_stamp; returns the current remote time stamp as a
            % double.
            r = this.send('t');
            ts = str2double(r);
        end
        
% 
%       function result = record(this)
%           import py.msgpack.loads
%           import org.zeromq.ZMQ;

%           sub_port=this.send('SUB_PORT');
%           sub = this.context.socket(ZMQ.SUB);
%           uri=sprintf('tcp://%s:%s', this.hostname, sub_port);
%           sub.connect(uri);
%           sub.subscribe('pupil.');
%           
%           for i=1:10
%               msg=sub.recvStr
%               data=sub.recvStr
%               c=char(data);
%               size(c)
%               loads(c)
%              % udata=msgpack('unpacker',data)

%           end
%           sub.close();
%           
%           result=data;
%       end
    end
end
