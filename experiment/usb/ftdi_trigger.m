classdef ftdi_trigger < handle
    %FTDI_TRIGGER - Use an FT232RQ USB-Serial converter's RTS pin to generate
    %TTL pulses.
    % t = ftdi_trigger('PORT') contstructs a ftdi_trigger object
    % PORT can be '/dev/tty.usbserial-FTA6FYZ4' (MacOS)
    %             'COM3' (Windows)
    % The trigger output is the ~RTS line of a FTDI232RQ USB-Serial
    % converter.
    %   
    % SEE ALSO: SET, TOGGLE, PULSE, STATUS
    properties (Access = protected)
        port
        current_value
    end
   
    methods
        function this = ftdi_trigger(portname)
            this.port = serial(portname, 'BaudRate', 4800);
            fopen(this.port);
            this.set(0);
        end
        
        function delete(this)
            fclose(this.port);
            delete(this.port);
        end
        
        function set(this,value)
            % SET - set output level
            % set(level) sets the output level to a logical 0 or 1
            assert(numel(value)==1);
            value=logical(value);
            pin = 'RequestToSend';
            if value   % the RTS output pin is inverted
                level='off';
            else
                level='on';
            end
            set(this.port,pin,level);
            this.current_value=value;
        end
        
        function toggle(this)
            %TOGGLE - toggle the output level
            this.set(~this.current_value);
        end
        
        function pulse(this,duration_seconds)
            %PULSE - generate a pulse
            %pulse(duration_in_seconds) toggles the output, waits for the
            %specified duration to elapse, then toggles back. 
            %Depending on the initial output level positive or negative
            %pulses can be generated.
            this.toggle();
            pause(duration_seconds);
            this.toggle();
        end
        
        function output_value = status(this)
            %STATUS - returns the current logical level of the output
            output_value = this.current_value;
        end
    end
end



