classdef (Abstract) lsldert_abstract_client < handle
    
    methods (Abstract)
        [result,telapsed]=send(this, msg)
    end
    
    methods
        function [result,telapsed] = beep(this, freq, duration, marker)
            % BEEP - play a beep on the raspberry
            % [result,telapsed] = obj.beep(freq, duration);
            % [result,telapsed] = obj.beep(880, 0.6); % plays a 880 Hz tone for 0.6
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
            [result,telapsed] = this.send(str);
        end
        
        function [result,telapsed] = pulseIR(this, mask, duration, marker)
            % PULSEIR - pulse one or two IR LEDs for fNIRS synchronization
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
            [result,telapsed] = this.send(str);
        end
        
        function [result,telapsed] = digitalout(this, value, marker)
            % DIGITALOUT - set digital output to high or low level
            if nargin < 3
                marker = sprintf("DOUT_%d", value);
            end
            str = sprintf("D %d %s", value, marker);
            [result,telapsed] = this.send(str);
        end
        
        function [result,telapsed] = set_digitalin_marker(this, marker)
            % SET_DIGITALIN_MARKER - set the marker name that is sent on digital in events
            str = sprintf("M %s", marker);
            [result,telapsed] = this.send(str);
        end
    end
    
    methods (Static)
        function zmq_write_multi(socket, cmd, varargin)
            import org.zeromq.*;
            
            % send command string followed by any data frames in varagin
            msg=ZMsg;
            msg.addString(sprintf('%s\0',cmd));
            ii=1;
            while ii <= numel(varargin)
                % transform the input data into a row-major
                % (lexicographically) ordered vector. This is the common
                % ordering for programming languages like C, C++,
                % Python (NumPy) and Pascal.
                if ischar(varargin{ii}) || isstring(varargin{ii})
                    frame=ZFrame(sprintf('%s\0',varargin{ii}));
                else
                    assert(isnumeric(varargin{ii}),'zmq_write_multi: expect numeric data');
                    data=reshape(transpose(varargin{ii}),1,[]);
                    bytes=typecast(data,'uint8');
                    frame=ZFrame(bytes);
                end
                msg.add(frame);
                ii=ii+1;
            end
            msg.send(socket);
        end
    end
end